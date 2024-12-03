//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MOIRE_POTENTIAL_H
#define QMCPLUSPLUS_MOIRE_POTENTIAL_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class MoirePotential : public OperatorBase
{
public:
  MoirePotential(ParticleSet& P)
  {
    setEnergyDomain(POTENTIAL);
    oneBodyQuantumDomain(P);
  };
  ~MoirePotential(){};
  Return_t evaluate(ParticleSet& P) override;
  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& P, TrialWaveFunction& psi) override;
  // ---- begin required overrides
  void resetTargetParticleSet(ParticleSet& P) override {APP_ABORT("not implemented");};
  std::string getClassName() const override {return "moire";};
  // required overrides end ----
private:
  RealType amoire, vmoire, phi, theta;
  std::vector<TinyVector<RealType, 3>> gvecs;
};
} // qmcplusplus
#endif
