#ifndef QMCPLUSPLUS_SPIN_ORIENTATION_H
#define QMCPLUSPLUS_SPIN_ORIENTATION_H

#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class SpinOrientation : public OperatorBase
{
public:
  SpinOrientation(const ParticleSet& ions, ParticleSet& P);

  std::string getClassName() const override { return "SpinOrientation"; }
  bool put(xmlNodePtr cur) override {         // read input xml node, required
    // add this estimator to stat.h5
    update_mode_.set(COLLECTABLE, 1);
    return true;
  }
  bool get(std::ostream& os) const override // class description, required
  {
    os << "SpinOrientation: " << name_;
    return true;
  }

  Return_t evaluate(ParticleSet& P) override;

  // record h5_index
  void registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const override;
  // allocate multiple columns
  void addObservables(PropertySetType& plist, BufferType& collectables) override;

  //should be empty for Collectables interface
  void resetTargetParticleSet(ParticleSet& P) override {}
  void setObservables(PropertySetType& plist) override {}
  void setParticlePropertyList(PropertySetType& plist, int offset) override {}
  // pure virtual functions require overrider
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) final { return std::make_unique<SpinOrientation>(*this); }

private:
  ParticleSet& tpset; // reference to target particle set
  const size_t npart;
  const size_t natom;
  const int itab;
  int h5_index; // index of this estimator in the collectables carried by target pset

}; // SpinOrientation

} // namespace qmcplusplus
#endif
