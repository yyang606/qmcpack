#include "Message/Communicate.h"
#include "QMCApp/QMCAppBase.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "LongRange/LRHandlerBase.h"
#include "LongRange/LPQHIBasis.h"
#include "LongRange/LRBreakup.h"
//#include <limits>

using namespace qmcplusplus;
using namespace std;

typedef QMCTraits::RealType RealType;


RealType fk(RealType k)
{ // Fourier representation of radial potential to break up
  return 4*M_PI/pow(k, 2);
}


RealType evaluateFklr(RealType k, vector<RealType> coefs, LPQHIBasis basis)
{
  RealType rc = basis.get_rc();
  RealType volume = basis.get_CellVolume();
  RealType fkval = fk(k);  //Natoli-Ceperley version //*cos(k*rc)/volume;
  for (int n=0; n<basis.NumBasisElem(); n++)
  {
    fkval -= coefs[n]*basis.c(n, k);
  }
  return fkval;
}


RealType sk0(RealType k, RealType kf)
{
  // handle boundary cases
  if (k<numeric_limits<RealType>::min()) return 0.0;
  if (k>=2*kf) return 1.0;
  RealType sk0_val = 3./4*(k/kf) - 1./16*pow(k/kf, 3);
  return sk0_val;
}


RealType sk(RealType k, RealType rs, RealType kf)
{
  RealType rho = 3./(4*M_PI*pow(rs, 3));
  RealType sk0_val = sk0(k, kf);
  RealType ak = 2*rho*pow(sk0_val, 2)*(8*M_PI/pow(k, 4));
  RealType sk_val = sk0_val*pow(1.+ak, -0.5);
  return sk_val;
}


RealType uk(RealType k, RealType rs, RealType kf)
{
  RealType rho = 3./(4*M_PI*pow(rs, 3));
  RealType sk0_val = sk0(k, kf);
  RealType sk_val = sk(k, rs, kf);
  RealType uk_val = (1./sk_val-1./sk0_val)/(2*rho);
  return uk_val;
}


RealType evaluateUklr(RealType k, RealType rs, RealType kf
  , vector<RealType> coefs, LPQHIBasis basis)
{
  RealType rc = basis.get_rc();
  RealType volume = basis.get_CellVolume();
  RealType val = uk(k, rs, kf);
  for (int n=0; n<basis.NumBasisElem(); n++)
  {
    val += coefs[n]*basis.c(n, k);
  }
  return val;
}


xmlNodePtr find(const char* expression, xmlXPathContextPtr context)
{ // find the first node matching xpath expression in the given context
  OhmmsXPathObject xpath(expression, context);
  if (xpath.size() != 1)
  {
    APP_ABORT("expected 1 " << expression << " found " << xpath.size());
  }
  return xpath[0];
}


int main(int argc, char **argv)
{
  OHMMS::Controller->initialize(argc, argv);
  // kgrid setup:
  //  0 < k < kcut: exact kshell degeneracy
  //  kcut < k < kmax: continuum approximate kshell degeneracy
  //
  // potential setup:
  //  0 < r < rc: frsr is non-zero
  //  0 < k < kc: fklr is non-zero
  //
  // to keep this script simple, everything is hard-coded for 3D
  // it should be easy to learn and adapt to 2D
 
  // input electron density using Wigner-Seitz radius
  RealType rs = 1.25;
  // !!!! hard-code unpolarized gas Fermi kvector
  RealType kf = pow(9*M_PI/4, 1./3)/rs;
  app_log() << "rs, kf " << rs << " " << kf << endl;

  // hard code some default parameters
  //  (stolen from LongRange/LRHandlerTemp.h)
  RealType kcut_prefactor = 60*M_PI;  // kcut = prefactor*(1/L)
  RealType kmax_rc = 6000.;
  // use default kcrc from box //RealType kc_rc = 15.;
  // approximate k-point degeneracies for kcut < k < kmax
  int nknot = 15;  // number of spline knots used in local function
  int nk = 1024;   // number of points on linear grid for output
  
  if (argc<2) APP_ABORT("usage: " << argv[0] << " axes.xml");

  // step 1: read <simulationcell> from input
  Libxml2Document fxml;
  fxml.parse(argv[1]);
  xmlXPathContextPtr doc = fxml.getXPathContext();
  xmlNodePtr sc_node = find("//simulationcell", doc);

  // step 2: construct simulation box
  Uniform3DGridLayout box;
  LatticeParser parser(box);
  parser.put(sc_node);
  box.SetLRCutoffs();  // box.LR_rc, box.LR_kc
  RealType volume = box.Volume;
  app_log() << "  lattice volume = " << volume << endl;

  // step 3: set short-range long-range (LR) breakup parameters
  RealType rc, kc, kcut, kmax;
  rc = box.LR_rc;
  kc = box.LR_kc;
  kcut = kcut_prefactor*pow(volume, -1./3.);
  kmax = kmax_rc/rc;
  app_log() << " rc, kc, kcut, kmax = " << rc << " " << kc 
    << " " << kcut << " " << kmax << endl;

  // step 4: define short-range basis functions (rc, nknot)
  LPQHIBasis basis(box);
  basis.set_Lattice(box);  // allow basis to access volume
  basis.set_NumKnots(nknot);
  basis.set_rc(rc);
  app_log() << "  basis cell volume = " << basis.get_CellVolume() << endl;
  
  // step 5: define long-range basis functions (plane waves on KList)
  LRBreakup<LPQHIBasis> handler(basis);
  int max_kshell = handler.SetupKVecs(kc, kcut, kmax);  // handler.KList
  app_log() << "  max_kshell = " << max_kshell << endl;
  // handler.KList contains two columns: |k| and degeneracy

  // put target function on kgrid, perform optimized breakup
  vector<RealType> xk, coefs;
  coefs.resize(basis.NumBasisElem());
  xk.resize(handler.KList.size());  // xk is FT[-V(r)] i.e. -f(k)
  for (int ik=0; ik<handler.KList.size(); ik++)
  {
    RealType kmag = handler.KList[ik][0];
    xk[ik] = fk(kmag);  // Natoli-Ceperley version //*(-cos(kmag*rc)/volume);
    // Q/ why *-cos/vol. ?
    // A/ Kenneth P. Esler 2006 thesis (5.50)
    //xk[ik] = uk(kmag, rs, kf);
  }
  
  // solve for expansion coefficients
  cout.precision(40);
  RealType chisq = handler.DoBreakup(xk.data(), coefs.data());
  app_log() << "  LR breakup chi^2 = " << chisq << endl;

  // output long-range potential Vlr(k) from kmin to kc with spacing dk
  RealType dk = kc/nk;
  //  set kmin to smallest reciprocal lattice vector
  RealType kmin = numeric_limits<RealType>::max();
  for (int i=0; i<3; i++)
  {
    RealType b = 2*M_PI*sqrt( dot(box.Gv[i], box.Gv[i]) );
    if (b<kmin) kmin=b;
  }
  app_log() << "#VK_START#" << endl;
  for (int ik=0; ik<nk; ik++)
  {
    RealType kmag = kmin+ik*dk;
    app_log() << kmag << " " << evaluateFklr(kmag, coefs, basis) << endl;
    //app_log() << kmag << " " << evaluateUklr(kmag, rs, kf, coefs, basis) << endl;
  }
  app_log() << "#VK_STOP#" << endl;

  app_log() << "#SK_START#" << endl;
  for (int ik=0; ik<nk; ik++)
  {
    RealType kmag = kmin+ik*dk;
    app_log() << kmag << " " << sk(kmag, rs, kf) << endl;
  }
  app_log() << "#SK_STOP#" << endl;

  app_log() << "#UK_START#" << endl;
  for (int ik=0; ik<nk; ik++)
  {
    RealType kmag = kmin+ik*dk;
    app_log() << kmag << " " << uk(kmag, rs, kf) << endl;
  }
  app_log() << "#UK_STOP#" << endl;
  
  OHMMS::Controller->finalize();
  return 0;
}
