#include "MoirePotential.h"
//#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{
MoirePotential::MoirePotential(ParticleSet& elec, ParticleSet& ions)
  : targetPtcl(elec),
    sourcePtcl(ions),
    ideitab(elec.addTable(ions))
{
}

MoirePotential::Return_t MoirePotential::evaluate(ParticleSet& P)
{
  return 0.0;
}

OperatorBase* MoirePotential::makeClone(ParticleSet& P, TrialWaveFunction& psi)
{
  return new MoirePotential(*this);
}

} // qmcplusplus
