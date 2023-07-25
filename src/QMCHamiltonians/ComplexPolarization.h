//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_COMPLEXPOLARIZATION_H
#define QMCPLUSPLUS_COMPLEXPOLARIZATION_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class ComplexPolarization : public OperatorBase
{
public:
  ComplexPolarization(ParticleSet& P);

  bool put(xmlNodePtr cur) override final; // read input xml node, required

  Return_t evaluate(ParticleSet& P) override final;
  // enable forward walking
  void setParticlePropertyList(PropertySetType& plist, int offset) override final;

  // allocate multiple columns in scalar.dat
  void addObservables(PropertySetType& plist, BufferType& collectables) override final;
  // fill multiple columns in scalar.dat
  void setObservables(PropertySetType& plist) override final;

  // ---- begin required overrides
  // pure virtual functions require overrider
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override final;
  void resetTargetParticleSet(ParticleSet& P) override final {throw std::runtime_error("not implemented");};
  std::string getClassName() const override final {return "ComplexPolarization";};
  bool get(std::ostream& os) const override final; // class description, required
  // required overrides end ----

private:
  ParticleSet& tpset; // reference to target particle set
  const ParticleSet::ParticleLayout& lattice; // used for frac. coord.
  const size_t ndim;
  std::vector<RealType> values;
  //  myIndex: the index of this estimator in the property list in target pset

}; // ComplexPolarization

} // namespace qmcplusplus
#endif
