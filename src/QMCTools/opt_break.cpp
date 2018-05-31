#include "Message/Communicate.h"
#include "QMCApp/QMCAppBase.h"
#include "QMCApp/ParticleSetPool.h"
#include "LongRange/LRHandlerBase.h"
#include "LongRange/LPQHIBasis.h"
#include "LongRange/LRBreakup.h"

using namespace qmcplusplus;
using namespace std;

typedef QMCTraits::RealType RealType;


xmlNodePtr find(const char* expression, xmlXPathContextPtr context)
{ // find the first node matching xpath expression in the given context
  OhmmsXPathObject xpath(expression, context);
  if (xpath.size() != 1)
  {
    APP_ABORT("expected 1 " << expression << " found " << xpath.size());
  }
  return xpath[0];
}


RealType fk(RealType k)
{ // Fourier representation of radial potential to break up
  return 4*M_PI/pow(k, 2);
}


RealType evaluateFklr(RealType k, vector<RealType> coefs, LPQHIBasis basis)
{
  RealType rc = basis.get_rc();
  RealType volume = basis.get_CellVolume();
  RealType fkval = fk(k)*cos(k*rc)/volume;
  for (int n=0; n<basis.NumBasisElem(); n++)
  {
    fkval += coefs[n]*basis.c(n, k);
  }
  return fkval;
}


int main(int argc, char **argv)
{
  OHMMS::Controller->initialize(argc, argv);
  
  if (argc<2) APP_ABORT("usage: " << argv[0] << " qmc.in.xml");

  // step 1: read input, must contain <simulationcell> and <particleset name="e">
  Libxml2Document fxml;
  fxml.parse(argv[1]);
  xmlXPathContextPtr doc = fxml.getXPathContext();

  // find simulationcell
  xmlNodePtr sc_node = find("//simulationcell", doc);

  // find quantum ("e") particleset
  xmlNodePtr epset_node = find("//particleset[@name='e']", doc);

  // step 2: construct ptclPool::ParticleSet with SK
  ParticleSetPool ptclPool(NULL);
  ptclPool.putLattice(sc_node); // <simulationcell>
  ptclPool.put(epset_node); // <particleset>

  // step 3: construct short-range basis functions
  ParticleSet* P = ptclPool.getParticleSet("e");
  LPQHIBasis basis(P->LRBox);
  basis.set_Lattice(P->LRBox);
  basis.set_NumKnots(15);
  basis.set_rc(P->LRBox.LR_rc);
  RealType volume = basis.get_CellVolume();
  app_log() << "cell volume = " << volume << endl;
  
  // step 4: construct long-range Coulomb handler
  LRBreakup<LPQHIBasis> handler(basis);
  RealType rc, kc, kcut, kmax;
  rc = P->LRBox.LR_rc;
  kc = 15./rc;
  kcut = 60*M_PI*pow(volume, -1./3.);
  kmax = 6000./rc;
  app_log() << " rc, kc, kcut, kmax = " << rc << " " << kc 
    << " " << kcut << " " << kmax << endl;
  int max_kshell = handler.SetupKVecs(kc, kcut, kmax);  // handler.KList
  app_log() << "max_kshell = " << max_kshell << endl;
  // handler.KList contains two columns: |k| and degeneracy

  // put function on kgrid, perform optimized breakup
  vector<RealType> xk, coefs;
  xk.resize(handler.KList.size());  // xk is FT[-V(r)] i.e. -f(k)
  for (int ik=0; ik<handler.KList.size(); ik++)
  {
    RealType kmag = handler.KList[ik][0];
    xk[ik] = -fk(kmag)*cos(kmag*rc)/volume;  // why *cos/vol. ?
  }
  coefs.resize(basis.NumBasisElem());
  
  RealType chisq = handler.DoBreakup(xk.data(), coefs.data());
  app_log() << "  LR breakup chi^2 = " << chisq << endl;

  // output long-range potential Vlr(k)
  app_log() << "#VK_START#" << endl;
  KContainer KList = P->SK->KLists;
  for (int ks=0, ik=0; ks<max_kshell; ks++)
  {
    RealType kmag = sqrt(KList.ksq[ik]);
    app_log() << kmag << " "
      << evaluateFklr(kmag, coefs, basis) << endl;
    while (ik<KList.kshell[ks+1] && KList.kpts_cart.size()) ik++;
  }
  app_log() << "#VK_STOP#" << endl;
  
  OHMMS::Controller->finalize();
  return 0;
}
