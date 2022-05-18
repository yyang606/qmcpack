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
  int npw;
  OhmmsAttributeSet attrib;
  attrib.add(npw, "size");
  attrib.put(cur);

  // extract npw k-points from container
  // kpts_cart is sorted by magnitude
  std::vector<PosType> kpts;
  kpts.resize(npw);
  kpts[0] = {0.0, 0.0, 0.0};
  for (int ik=1;ik<npw;ik++)
  {
    kpts[ik] = targetPtcl.SK->KLists.kpts_cart[ik-1];
  }
  PWOrbital* sposet = new PWOrbital(kpts);
  sposet->report("  ");
  return sposet;
}

} // qmcplusplus
