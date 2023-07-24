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

#ifndef QMCPLUSPLUS_VECTORPAIRCORR_H
#define QMCPLUSPLUS_VECTORPAIRCORR_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class VectorPairCorr : public OperatorBase
{
public:
  VectorPairCorr(ParticleSet& P);

  bool put(xmlNodePtr cur) override final; // read input xml node, required

  Return_t evaluate(ParticleSet& P) override final;

  // allocate multiple columns in scalar.dat (plist) or stat.h5 (collectables)
  void addObservables(PropertySetType& plist, BufferType& collectables) override final;
  // fill columns in scalar.dat
  void setObservables(PropertySetType& plist) override final {};
  // fill datasets in stat.h5
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const override final;

  // ---- begin required overrides
  // pure virtual functions require overrider
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override final;
  void resetTargetParticleSet(ParticleSet& P) override final {APP_ABORT("not implemented");};
  std::string getClassName() const override final {return "VectorPairCorr";};
  bool get(std::ostream& os) const override final; // class description, required
  // required overrides end ----

private:
  ParticleSet& tpset; // reference to target particle set
  const ParticleSet::ParticleLayout& lattice; // used for frac. coord.
  const size_t ndim;
  //  my_index_: the index of this estimator in the property list in target pset
  const int d_aa_ID_;
  TinyVector<int, DIM> grid;
  TinyVector<int, DIM> gdims;
  int npoints;

}; // VectorPairCorr

} // namespace qmcplusplus
#endif
