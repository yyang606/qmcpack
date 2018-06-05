#include "Message/Communicate.h"
#include "QMCApp/QMCAppBase.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "ParticleIO/ParticleLayoutIO.h"
#include "LongRange/LRHandlerBase.h"
#include "LongRange/LPQHIBasis.h"
#include "LongRange/LRBreakup.h"
#include "LongRange/KContainer.h"
//#include <limits>

using namespace qmcplusplus;
using namespace std;

typedef QMCTraits::RealType RealType;


// LR break Coulomb potential
RealType eval_vk(RealType k)
{ // bare Coulomb potential
  return 4*M_PI/pow(k, 2);
}


RealType eval_vklr(RealType k, vector<RealType> coefs, LPQHIBasis basis)
{ // Natoli long-range part of Coulomb potential
  RealType val = eval_vk(k);
  for (int n=0; n<basis.NumBasisElem(); n++)
  {
    val -= coefs[n]*basis.c(n, k);
  }
  return val;
}


// LR break Jastrow potential: need to set up Gaskell U(k) first
RealType eval_sk0(RealType k, RealType kf)
{ // wf = ideal Fermi gas determinant
  // handle boundary cases
  if (k<numeric_limits<RealType>::min()) return 0.0;
  if (k>=2*kf) return 1.0;
  RealType sk0 = 3./4*(k/kf) - 1./16*pow(k/kf, 3);
  return sk0;
}


RealType eval_sk(RealType k, RealType rs, RealType kf)
{ // wf = ideal Fermi gas determinant * Gaskell RPA Jastrow
  RealType rho = 3./(4*M_PI*pow(rs, 3));
  RealType sk0 = eval_sk0(k, kf);
  RealType ak = 2*rho*pow(sk0, 2)*(8*M_PI/pow(k, 4));
  RealType sk_val = sk0*pow(1.+ak, -0.5);
  return sk_val;
}


RealType eval_uk(RealType k, RealType rs, RealType kf)
{ // Gaskell RPA Jastrow potential
  RealType rho = 3./(4*M_PI*pow(rs, 3));
  RealType sk0 = eval_sk0(k, kf);
  RealType sk = eval_sk(k, rs, kf);
  RealType uk = (1./sk-1./sk0)/(2*rho);
  return uk;
}


RealType eval_uklr(RealType k, RealType rs, RealType kf
  , vector<RealType> coefs, LPQHIBasis basis)
{ // Natoli long-range part of RPA Jastrow potential
  RealType val = eval_uk(k, rs, kf);
  for (int n=0; n<basis.NumBasisElem(); n++)
  {
    val -= coefs[n]*basis.c(n, k);
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
  if (argc<3) APP_ABORT("usage: " << argv[0] << " axes.xml rs");
 
  // input electron density using Wigner-Seitz radius
  RealType rs = stof(argv[2]);
  // !!!! hard-code unpolarized gas Fermi kvector
  RealType kf = pow(9*M_PI/4, 1./3)/rs;
  app_log() << "rs, kf " << rs << " " << kf << endl;
  // electron density in Hartree atomic units
  RealType rho = 3./(4*M_PI*pow(rs, 3));

  // hard code some default parameters
  //  (stolen from LongRange/LRHandlerTemp.h)
  RealType kcut_prefactor = 60*M_PI;  // kcut = prefactor*(1/L)
  RealType kmax_rc = 6000.;
  // use default kcrc from box //RealType kc_rc = 15.;
  // approximate k-point degeneracies for kcut < k < kmax
  int nknot = 15;  // number of spline knots used in local function
  int nk = 1024;   // number of points on linear grid for output

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

  // step 6: put target function on kgrid
  vector<RealType> vkvals, ukvals;
  vkvals.resize(handler.KList.size());
  ukvals.resize(handler.KList.size());
  for (int ik=0; ik<handler.KList.size(); ik++)
  {
    RealType kmag = handler.KList[ik][0];
    vkvals[ik] = eval_vk(kmag);
    ukvals[ik] = eval_uk(kmag, rs, kf);
  }

  // step 7: perform optimized breakup i.e. solve for expansion coefficients
  vector<RealType> vkcoefs, ukcoefs;
  vkcoefs.resize(basis.NumBasisElem());
  ukcoefs.resize(basis.NumBasisElem());
  RealType vkchisq = handler.DoBreakup(vkvals.data(), vkcoefs.data());
  RealType ukchisq = handler.DoBreakup(ukvals.data(), ukcoefs.data());
  cout.precision(40);
  app_log() << "  Vk LR breakup chi^2 = " << vkchisq << endl;
  app_log() << "  Uk LR breakup chi^2 = " << ukchisq << endl;

  RealType dk = kc/nk;
  RealType kmin = 1e-1;

  // output Coulomb Vklr for debugging
  app_log() << "#VK_START#" << endl;
  for (int ik=0; ik<nk; ik++)
  {
    RealType kmag = kmin+ik*dk;
    app_log() << kmag << " " << eval_vklr(kmag, vkcoefs, basis) << endl;
  }
  app_log() << "#VK_STOP#" << endl;

  // step 8: perform Holzmann 2016 finite size correction (FSC)
  //  use eq. (30) for potential and eq. (35) for kinetic

  // output integrand for the long-range part of potential and kinetic
  //  use fine kgrid for 1D quadrature
  app_log() << "#VFSC_START#" << endl;
  for (int ik=0; ik<nk; ik++)
  {
    RealType kmag = kmin+ik*dk;
    RealType vklr, sk, integrand;
    vklr = eval_vklr(kmag, vkcoefs, basis);
    sk = eval_sk(kmag, rs, kf);
    integrand = pow(kmag, 2)/(2*M_PI*M_PI)* 0.5*vklr*sk;  // isotropic 3D -> 1D
    app_log() << kmag << " " << integrand << endl;
  }
  app_log() << "#VFSC_STOP#" << endl;

  app_log() << "#TFSC_START#" << endl;
  for (int ik=0; ik<nk; ik++)
  {
    RealType kmag = kmin+ik*dk;
    RealType uk, uklr, sk, integrand, ukpiece;
    uk = eval_uk(kmag, rs, kf);
    sk = eval_sk(kmag, rs, kf);
    uklr = eval_uklr(kmag, rs, kf, ukcoefs, basis);
    //integrand = pow(kmag, 2)/(2*M_PI*M_PI)* 0.5*rho*pow(kmag, 2)*uklr*(2*uk-uklr)*sk;
    // group k^2U(k) together to avoid numerical instability
    ukpiece = (pow(kmag, 2)*uklr) * ((2*pow(kmag, 2)*uk-pow(kmag, 2)*uklr));
    integrand = 0.5/(2*M_PI*M_PI)*rho*ukpiece*sk;

    app_log() << kmag << " " << integrand << endl;
  }
  app_log() << "#TFSC_STOP#" << endl;

  // output long-range part of potential and kinetic sums
  //  use only kvectors included in the simulation
  KContainer kvecs;
  kvecs.UpdateKLists(box, kc);

  RealType vsum = 0.0;
  RealType tsum = 0.0;
  for (int ik=0; ik<kvecs.ksq.size(); ik++)
  {
    RealType kmag, vklr, sk, uk, uklr;
    kmag = sqrt(kvecs.ksq[ik]);

    vklr = eval_vklr(kmag, vkcoefs, basis);
    sk   = eval_sk(kmag, rs, kf);
    vsum += 0.5*vklr*sk;

    uk   = eval_uk(kmag, rs, kf);
    uklr = eval_uklr(kmag, rs, kf, ukcoefs, basis);
    tsum += 0.5*rho*pow(kmag, 2)*uklr*(2*uk-uklr)*sk;
  }
  app_log() << "  vsum = " << vsum/volume << endl;
  app_log() << "  tsum = " << tsum/volume << endl;
  
  OHMMS::Controller->finalize();
  return 0;
}
