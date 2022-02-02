#ifndef QMCPLUSPLUS_GAUSSIAN_ORBITAL_BUILDER_H
#define QMCPLUSPLUS_GAUSSIAN_ORBITAL_BUILDER_H

#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "GaussianOrbitalSet.h"

namespace qmcplusplus
{
class GaussianOrbitalBuilder : public SPOSetBuilder
{
public:
  typedef std::map<std::string, ParticleSet*> PtclPoolType;

  GaussianOrbitalBuilder(ParticleSet& els, PtclPoolType& psets, Communicate* comm, xmlNodePtr cur);
  ~GaussianOrbitalBuilder();

  SPOSet* createSPOSetFromXML(xmlNodePtr cur);
private:
  ParticleSet& targetPtcl;
  std::string sourceName;
  PtclPoolType& particleSets;
};
} // qmcplusplus
#endif
