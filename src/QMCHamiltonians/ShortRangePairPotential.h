#ifndef QMCPLUSPLUS_SHORTRANGE_PAIR_POTENTIAL_H
#define QMCPLUSPLUS_SHORTRANGE_PAIR_POTENTIAL_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
struct ShortRangePairPotential : public OperatorBase
{
  //data members
  RealType amplitute, sigma;

  //construction/destruction
  ShortRangePairPotential(ParticleSet& P)
  {
    set_energy_domain(potential);
    two_body_quantum_domain(P);
  }

  ~ShortRangePairPotential() {}

  //unneeded interface functions
  void resetTargetParticleSet(ParticleSet& P) {}

  //standard interface functions
  bool put(xmlNodePtr cur);
  bool get(std::ostream& os) const;
  OperatorBase* makeClone(ParticleSet& P, TrialWaveFunction& psi);

  //functions for physical (hamiltonian component) estimator
  Return_t evaluate(ParticleSet& P);
  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy) { return evaluate(P); }
};
} // namespace qmcplusplus
#endif
