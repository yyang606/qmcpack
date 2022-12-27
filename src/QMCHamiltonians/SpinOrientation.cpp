#include "SpinOrientation.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

SpinOrientation::SpinOrientation(const ParticleSet& ions, ParticleSet& P) :
  tpset(P),
  npart(P.getTotalNum()),
  natom(ions.getTotalNum()),
  itab(P.addTable(ions))
{
};

SpinOrientation::Return_t SpinOrientation::evaluate(ParticleSet& P)
{
  const auto& d_table = P.getDistTableAB(itab);
  RealType wgt = t_walker_->Weight;
  for (int i = 0; i < npart; i++)
  {
    const auto& dist = d_table.getDistRow(i);
    const int iat = std::min_element(dist.begin(), dist.begin()+natom) - dist.begin();
    P.Collectables[h5_index + iat] += P.spins[iat] * wgt;
  }

  value_ = 0.0; // Value is no longer used
  return value_;
}

// ----------------------------- boiler plate -----------------------------
void SpinOrientation::addObservables(PropertySetType& plist, BufferType& collectables)
{
  // make room in h5 file and save its index
  h5_index = collectables.size();
  std::vector<RealType> tmp(natom);
  collectables.add(tmp.begin(), tmp.end());
}

void SpinOrientation::registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  std::vector<int> ndim(1, natom);
  h5desc.emplace_back(hdf_path{name_});
  auto& h5o = h5desc.back();
  h5o.set_dimensions(ndim, h5_index);
}

} // namespace qmcplusplus
