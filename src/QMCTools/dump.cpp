#include "Message/Communicate.h"
#include "QMCApp/QMCAppBase.h"
#include "QMCApp/ParticleSetPool.h"
#include "LongRange/LRCoulombSingleton.h"

using namespace qmcplusplus;
using namespace std;

typedef QMCTraits::RealType RealType;


xmlNodePtr find(const char* expression, xmlXPathContextPtr context)
{ // find the first node matching xpath expression in the given context
  OhmmsXPathObject xpath(expression, context);
  if (xpath.size() != 1) APP_ABORT("expected 1 " << expression << " found " << xpath.size());
  return xpath[0];
}


int main(int argc, char **argv)
{
  OHMMS::Controller->initialize(argc,argv);
  
  if (argc<2) APP_ABORT("usage: ./dump qmc.in.xml");

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
  
  // step 3: construct long-range Coulomb handler
  LRCoulombSingleton::LRHandlerType* AA;
  ParticleSet* P = ptclPool.getParticleSet("e");
  AA = LRCoulombSingleton::getHandler(*P);
  RealType rcut = AA->get_rc();
  cout << "LR rcut = " << rcut << endl;
  cout << "LR kcut = " << AA->get_kc() << endl;

  // output long-range potential Vlr(k)
  KContainer klist = P->SK->KLists;

  app_log() << "#VK_START#" << endl;
  for (int ks=0; ks<klist.kshell.size()-1; ks++)
  {
    app_log() << sqrt(klist.ksq[klist.kshell[ks]]) << " "
      << AA->Fk_symm[ks] << endl;
  }
  app_log() << "#VK_STOP#" << endl;
  
  // step 4: construct short-range Colomb spline
  //LRCoulombSingleton::GridType* grid(0);
  LRCoulombSingleton::RadFunctorType* rVs;
  rVs = LRCoulombSingleton::createSpline4RbyVs(AA, rcut, 0);

  // output short-range potential r*Vsr(r)
  int npts = 1001;  // number of grid points
  RealType r0=0, r1=rcut;
  RealType dr = (r1-r0)/npts;  // lineargrid spacing

  app_log() << "#RVR_START#" << endl;
  for (int i=0; i<npts; i++)
  {
    RealType r = r0 + dr*i;
    app_log() << r << " " << rVs->splint(r) << endl;
  }
  app_log() << "#RVR_STOP#" << endl;

  OHMMS::Controller->finalize();
  return 0;
}
