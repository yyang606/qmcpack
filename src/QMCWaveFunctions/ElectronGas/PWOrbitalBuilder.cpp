#include "PWOrbitalBuilder.h"
#include "OhmmsData/AttributeSet.h"
#include "LongRange/StructFact.h"
#include "PWOrbital.h"

namespace qmcplusplus
{

PWOrbitalBuilder::PWOrbitalBuilder(ParticleSet& els, Communicate* comm, xmlNodePtr cur)
  : SPOSetBuilder("PW", comm),
    targetPtcl(els)
{}

SPOSet* PWOrbitalBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  int npw, norb;
  PosType twist(0.0);
  OhmmsAttributeSet attrib;
  attrib.add(norb, "size");
  attrib.add(twist, "twist");
  attrib.put(cur);

  auto lattice = targetPtcl.Lattice;

#ifdef QMC_COMPLEX
  npw = norb;
#else
  npw = std::ceil((norb+1.0)/2);
#endif
  targetPtcl.setTwist(twist);
  PosType tvec = lattice.k_cart(twist);
  app_log() << "twist fraction = " << twist << std::endl;
  app_log() << "twist cartesian = " << tvec << std::endl;

  // extract npw k-points from container
  // kpts_cart is sorted by magnitude
  std::vector<PosType> kpts(npw);
  KContainer klists;
  RealType kcut = lattice.LR_kc; // to-do: reduce kcut to >~ kf
  klists.UpdateKLists(lattice, kcut, twist, lattice.ndim);

  kpts[0] = tvec;
  for (int ik=1;ik<npw;ik++)
  {
    kpts[ik] = tvec+klists.kpts_cart[ik-1];
  }
  PWOrbital* sposet = new PWOrbital(kpts);
  sposet->report("  ");
  return sposet;
}

} // qmcplusplus
