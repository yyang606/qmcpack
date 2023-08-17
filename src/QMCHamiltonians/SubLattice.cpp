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

#include "SubLattice.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/DistanceTable.h"
#include "asp.h"

namespace qmcplusplus
{
SubLattice::SubLattice(ParticleSet& P, ParticleSet& source) :
  tpset(P),
  spset(source),
  lattice(P.getLattice()),
  ndim(lattice.ndim),
  nelec(P.getTotalNum()),
  idtable(P.addTable(source))
{
  update_mode_.set(COLLECTABLE, 1);
  auto sspec = source.getSpeciesSet();
  const int nspec = sspec.size();
  firsts.resize(nspec);
  lasts.resize(nspec);
  for (int i=0;i<sspec.size();i++)
  {
    int ns = source.first(i);
    int n1 = source.last(i);
    firsts[i] = ns;
    lasts[i] = n1;
    //app_log() << i << " " << ns << " " << n1 << std::endl;
  }
  rij.resize(nelec);
  for (int i=0; i<nelec; ++i)
    rij[i].resize(nelec);
  temp_rij = create_matrix(nelec, nelec);
  ij_map.resize(nelec);
  for (int i=0; i<ij_map.size();i++)
    ij_map[i] = i;
};

bool SubLattice::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attrib;
  attrib.put(cur);
  get(app_log());
  return true;
}

bool SubLattice::get(std::ostream& os) const
{ // class description
  os << "SubLattice: " << name_ << std::endl;
  return true;
}

void SubLattice::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = collectables.current();
  std::vector<int> tmp(nelec);
  collectables.add(tmp.begin(), tmp.end());
}

void SubLattice::registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  // one scalar per lattice site (i.e. source particle)
  std::vector<int> ng(1);
  ng[0] = nelec;
  h5desc.emplace_back(hdf_path{name_});
  auto& oh = h5desc.back();
  oh.set_dimensions(ng, my_index_);
}



SubLattice::Return_t SubLattice::evaluate(ParticleSet& P)
{
  RealType wgt = t_walker_->Weight;

  const auto& d_table = P.getDistTableAB(idtable);
  // !!!! HACK: use all electrons
  const int first_tar = 0;
  const int last_tar = nelec;

  std::vector<RealType> dists_min(nelec);
  std::fill(dists_min.begin(), dists_min.end(), 999999999999.0);
  std::vector<int> min_ids(nelec);
  std::vector<int> min_specs(nelec);
  for (int ispec=0; ispec<firsts.size(); ispec++)
  {
    // extract distance table section
    const int first_src = firsts[ispec];
    const int last_src = lasts[ispec];
    for (int iel=first_tar;iel<last_tar;iel++)
    {
      const int i = iel-first_tar;
      const auto& dists = d_table.getDistRow(iel);
      for (int jat=first_src;jat<last_src;jat++)
      {
        const int j = jat-first_src;
        const RealType r = dists[jat];
        rij[i][j] = r;
        temp_rij[i][j] = r; // copy rij to be scrambled by asp
      }
    }
    // assign electron to sites
    for (int i=0; i<ij_map.size();i++)
      ij_map[i] = i;
    asp(ij_map.size(), temp_rij, ij_map.data());
    // update min
    for (int iel=first_tar;iel<last_tar;iel++)
    {
      const int i = iel-first_tar;
      const int j = ij_map[i];
      const auto r = rij[i][j];
      if (r < dists_min[i])
      {
        dists_min[i] = r;
        min_ids[i] = j+first_src;
        min_specs[i] = ispec;
      }
    }
  }
  for (int i=0;i<nelec;i++)
  {
    P.Collectables[my_index_+i] += wgt*min_specs[i];
  }
  return value_;
}

std::unique_ptr<OperatorBase> SubLattice::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<SubLattice>(*this);
}

} // namespace qmcplusplus
