#include "GaussianOrbitalBuilder.h"

namespace qmcplusplus
{

GaussianOrbitalBuilder::GaussianOrbitalBuilder(ParticleSet& els, PtclPoolType& psets, Communicate* comm, xmlNodePtr cur)
  : SPOSetBuilder("Gaussian", comm),
    targetPtcl(els)
{}

GaussianOrbitalBuilder::~GaussianOrbitalBuilder(){}

SPOSet* GaussianOrbitalBuilder::createSPOSetFromXML(xmlNodePtr cur)
{
  GaussianOrbitalSet* sposet = new GaussianOrbitalSet();
  return sposet;
}

} // qmcplusplus
