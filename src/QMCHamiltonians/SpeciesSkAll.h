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

#ifndef QMCPLUSPLUS_SPECIESSKALL_H
#define QMCPLUSPLUS_SPECIESSKALL_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class SpeciesSkAll : public OperatorBase
{
public:
  SpeciesSkAll(ParticleSet& P);

  bool put(xmlNodePtr cur) override final; // read input xml node, required

  Return_t evaluate(ParticleSet& P) override final;

  // allocate multiple columns in scalar.dat (plist) or stat.h5 (collectables)
  void addObservables(PropertySetType& plist, BufferType& collectables) override final;
  // fill columns in scalar.dat
  void setObservables(PropertySetType& plist) override final {}
  // fill datasets in stat.h5
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const override final;
  // !!!! Must override to avoid memory corruption
  void setParticlePropertyList(PropertySetType& plist, int offset) override final {}

  // ---- begin required overrides
  // pure virtual functions require overrider
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override final;
  void resetTargetParticleSet(ParticleSet& P) override final {APP_ABORT("not implemented");};
  std::string getClassName() const override final {return "SpeciesSkAll";};
  bool get(std::ostream& os) const override final; // class description, required
  // required overrides end ----

  int gen_pair_id(const int ig, const int jg, const int ns) const;

private:
  const ParticleSet::ParticleLayout& lattice; // used for frac. coord.
  const size_t ndim;
  //  my_index_: the index of this estimator in the property list in target pset
  const SpeciesSet& species;
  const size_t nspec, npair, npoints;
}; // SpeciesSkAll

} // namespace qmcplusplus
#endif
