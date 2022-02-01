#include "GaussianOrbitalBuilder.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

GaussianOrbitalBuilder::GaussianOrbitalBuilder(ParticleSet& els, PtclPoolType& psets, Communicate* comm, xmlNodePtr cur)
  : SPOSetBuilder("Gaussian", comm),
    targetPtcl(els),
    particleSets(psets)
{}

GaussianOrbitalBuilder::~GaussianOrbitalBuilder(){}

SPOSet* GaussianOrbitalBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  RealType cexpo;
  int ndim;
  std::string sourceName;
  OhmmsAttributeSet attrib;
  attrib.add(cexpo, "c");
  attrib.add(ndim, "ndim");
  attrib.add(sourceName, "source");
  attrib.put(cur);
  ParticleSet* sourcePtcl = particleSets[sourceName];
  GaussianOrbitalSet* sposet = new GaussianOrbitalSet(targetPtcl, *sourcePtcl, cexpo, ndim);
  sposet->report("  ");
  return sposet;
}

} // qmcplusplus
