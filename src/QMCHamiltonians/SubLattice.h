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

#ifndef QMCPLUSPLUS_SUBLATTICE_H
#define QMCPLUSPLUS_SUBLATTICE_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class SubLattice : public OperatorBase
{
public:
  SubLattice(ParticleSet& P, ParticleSet& source);

  bool put(xmlNodePtr cur) override final; // read input xml node, required

  Return_t evaluate(ParticleSet& P) override final;

  void addObservables(PropertySetType& plist, BufferType& collectables) override final;
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const override final;

  // ---- begin required overrides
  // pure virtual functions require overrider
  void setObservables(PropertySetType& plist) override final {}
  void setParticlePropertyList(PropertySetType& plist, int offset) override final {}
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override final;
  void resetTargetParticleSet(ParticleSet& P) override final {throw std::runtime_error("not implemented");};
  std::string getClassName() const override final {return "SubLattice";};
  bool get(std::ostream& os) const override final; // class description, required
  // required overrides end ----

private:
  ParticleSet& tpset; // reference to target particle set
  ParticleSet& spset; // reference to source particle set
  const ParticleSet::ParticleLayout& lattice;
  const size_t ndim;
  const size_t nelec;
  const int idtable;
  std::vector<int> firsts, lasts;
  std::vector<std::vector<RealType>> rij;
  RealType **temp_rij;  // interface to asp routine
  std::vector<long> ij_map;
}; // SubLattice

} // namespace qmcplusplus
#endif
