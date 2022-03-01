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
  RealType cexpo=0.1;
  int ndim=3;
  std::string sourceName="ion0";
  OhmmsAttributeSet attrib;
  attrib.add(cexpo, "c");
  attrib.add(ndim, "ndim");
  attrib.add(sourceName, "source");
  attrib.put(cur);
  ParticleSet* sourcePtcl = particleSets[sourceName];
  if (sourcePtcl == 0)
  {
    APP_ABORT("GaussianOrbital needs the source particleset");
  }
  GaussianOrbital* sposet = new GaussianOrbital(targetPtcl, *sourcePtcl, cexpo, ndim);
  sposet->report("  ");
  return sposet;
}

} // qmcplusplus
