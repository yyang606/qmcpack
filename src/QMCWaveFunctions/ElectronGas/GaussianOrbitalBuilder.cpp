#include "GaussianOrbitalBuilder.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

GaussianOrbitalBuilder::GaussianOrbitalBuilder(ParticleSet& els, PtclPoolType& psets, Communicate* comm, xmlNodePtr cur)
  : SPOSetBuilder("Gaussian", comm),
    targetPtcl(els)
{}

GaussianOrbitalBuilder::~GaussianOrbitalBuilder(){}

SPOSet* GaussianOrbitalBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  OhmmsAttributeSet attrib;
  attrib.add(cexpo, "c");
  attrib.put(cur);
  GaussianOrbitalSet* sposet = new GaussianOrbitalSet(cexpo);
  sposet->report("  ");
  return sposet;
}

} // qmcplusplus
