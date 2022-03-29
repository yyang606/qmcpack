#ifndef QMCPLUSPLUS_PW_ORBITAL_BUILDER_H
#define QMCPLUSPLUS_PW_ORBITAL_BUILDER_H

#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "PWOrbital.h"

namespace qmcplusplus
{
class PWOrbitalBuilder : public SPOSetBuilder
{
public:
  PWOrbitalBuilder(ParticleSet& els, Communicate* comm, xmlNodePtr cur);
  ~PWOrbitalBuilder(){}

  SPOSet* createSPOSetFromXML(xmlNodePtr cur);
private:
  ParticleSet& targetPtcl;
};
} // qmcplusplus
#endif
