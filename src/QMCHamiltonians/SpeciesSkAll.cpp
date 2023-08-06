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

#include "SpeciesSkAll.h"
#include "OhmmsData/AttributeSet.h"
#include "LongRange/StructFact.h"

namespace qmcplusplus
{
SpeciesSkAll::SpeciesSkAll(ParticleSet& P) :
  lattice(P.getLattice()),
  ndim(lattice.ndim),
  species(P.getSpeciesSet()),
  nspec(species.size()),
  npair((nspec*(nspec-1))/2+nspec),
  npoints(P.getSimulationCell().getKLists().numk),
  kvecs(P.getSimulationCell().getKLists().kpts_cart)
{
  update_mode_.set(COLLECTABLE, 1);
  bool periodic = lattice.SuperCellEnum != SUPERCELL_OPEN;
  if (not periodic)
  {
    std::ostringstream msg;
    msg << "NotImplementedError: SpeciesSkAll with open boundary\n";
    throw std::runtime_error(msg.str());
  }
};

bool SpeciesSkAll::put(xmlNodePtr cur)
{
  get(app_log());
  return true;
}

bool SpeciesSkAll::get(std::ostream& os) const
{ // class description
  os << "SpeciesSkAll: " << name_ << std::endl;
  os << "  ndim    = " << ndim << std::endl;
  os << "  npair   = " << npair << std::endl;
  return true;
}

void SpeciesSkAll::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = collectables.current();
  std::vector<RealType> tmp(npair*npoints+nspec*2*npoints);
  collectables.add(tmp.begin(), tmp.end());
}

void SpeciesSkAll::registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  // add lattice information
  Matrix<double> axes(ndim);
  for (int l=0;l<ndim;l++)
    for (int m=0;m<ndim;m++)
      axes(l,m) = lattice.R(l,m);
  file.push(hdf_path{name_});
  file.write(axes, "axes");
  file.pop();
  // Add k-point information
  h5desc.emplace_back(hdf_path{name_} / "kpoints");
  auto& obs_helper = h5desc.back();
  obs_helper.addProperty(kvecs, "value", file);
  // space for sofk
  std::vector<int> ng(2);
  ng[0] = npair;
  ng[1] = npoints;
  h5desc.emplace_back(hdf_path{name_}/"sofk");
  auto& oh = h5desc.back();
  oh.set_dimensions(ng, my_index_);
  // space for rhok
  ng[0] = nspec;
  ng[1] = 2*npoints; // real & imag
  h5desc.emplace_back(hdf_path{name_}/"rhok");
  auto& rhok_obs_helper = h5desc.back();
  rhok_obs_helper.set_dimensions(ng, my_index_+npair*npoints);
}

SpeciesSkAll::Return_t SpeciesSkAll::evaluate(ParticleSet& P)
{
  RealType wt = t_walker_->Weight;
 
  // sofk
  for (int ig=0;ig<nspec;ig++)
  {
    auto* restrict rki_re_ptr = P.getSK().rhok_r[ig];
    auto* restrict rki_im_ptr = P.getSK().rhok_i[ig];
    for (int jg=0;jg<nspec;jg++)
    {
      auto* restrict rkj_re_ptr = P.getSK().rhok_r[jg];
      auto* restrict rkj_im_ptr = P.getSK().rhok_i[jg];
      const int ipair = gen_pair_id(ig, jg, nspec);
      const int offset = my_index_+ipair*npoints;
      for (int k=0;k<npoints;k++)
      {
        P.Collectables[offset+k] += wt*(rki_re_ptr[k]*rkj_re_ptr[k] + rki_im_ptr[k]*rkj_im_ptr[k]);
      }
    }
  }
  // rhok
  for (int ig=0;ig<nspec;ig++)
  {
    auto* restrict rki_re_ptr = P.getSK().rhok_r[ig];
    auto* restrict rki_im_ptr = P.getSK().rhok_i[ig];
    const int rhok_offset = my_index_+npair*npoints+ig*2*npoints;
    for (int k=0;k<npoints;k++)
    {
      P.Collectables[rhok_offset+k] += wt*rki_re_ptr[k];
      P.Collectables[rhok_offset+npoints+k] += wt*rki_im_ptr[k];
    }
  }

  value_ = 0.0; // Value is no longer used in scalar.dat
  return value_;
}

int SpeciesSkAll::gen_pair_id(const int ig, const int jg, const int ns) const
{ // steal from PairCorr
  if (jg < ig)
    return ns * (ns - 1) / 2 - (ns - jg) * (ns - jg - 1) / 2 + ig;
  else
    return ns * (ns - 1) / 2 - (ns - ig) * (ns - ig - 1) / 2 + jg;
}

std::unique_ptr<OperatorBase> SpeciesSkAll::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<SpeciesSkAll>(*this);
}

} // namespace qmcplusplus