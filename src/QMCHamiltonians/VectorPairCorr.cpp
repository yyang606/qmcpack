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

#include "VectorPairCorr.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{
VectorPairCorr::VectorPairCorr(ParticleSet& P) :
  lattice(P.getLattice()),
  ndim(lattice.ndim),
  d_aa_ID_(P.addTable(P, DTModes::NEED_FULL_TABLE_ON_HOST_AFTER_DONEPBYP)),
  species(P.getSpeciesSet()),
  npair((species.size()*(species.size()-1))/2+species.size()),
  grid(-1)
{
  update_mode_.set(COLLECTABLE, 1);
  bool periodic = lattice.SuperCellEnum != SUPERCELL_OPEN;
  if (not periodic)
  {
    std::ostringstream msg;
    msg << "NotImplementedError: VectorPairCorr with open boundary\n";
    throw std::runtime_error(msg.str());
  }
  const size_t nspec=species.size();
  norms.resize(nspec, nspec);
  for (int s = 0; s < nspec; ++s)
    species_size.push_back(P.groupsize(s));
};

bool VectorPairCorr::put(xmlNodePtr cur)
{
  xmlNodePtr element = cur->xmlChildrenNode;
  while (element != NULL)
  {
    std::string ename((const char*)element->name);
    if (ename == "parameter")
    {
      const std::string name(getXMLAttributeValue(element, "name"));
      if (name == "grid")
      {
        putContent(grid, element);
      }
    }
    element = element->next;
  }
  // number of grid points
  npoints = 1;
  for (int l = 0; l < ndim; ++l)
  {
    if (grid[l] <= 0)
    {
      std::ostringstream msg;
      msg << "grid is required in VectorPairCorr\n";
      throw std::runtime_error(msg.str());
    }
    npoints *= grid[l];
  }
  gdims[0] = npoints / grid[0];
  for (int d = 1; d < ndim; ++d)
    gdims[d] = gdims[d - 1] / grid[d];
  // normalize
  for (int i=0;i<species.size();i++)
  {
    const size_t ni = species_size[i];
    for (int j=0;j<species.size();j++)
    {
      const size_t nj = species_size[j];
      const size_t npij =  i==j ? (ni*(nj-1))/2 : ni*nj;
      norms(i,j) = (RealType)npoints/npij; // 1/[average hit per bin]
    }
  }
  // report
  get(app_log());
  return true;
}

bool VectorPairCorr::get(std::ostream& os) const
{ // class description
  os << "VectorPairCorr: " << name_ << std::endl;
  os << "  ndim    = " << ndim << std::endl;
  os << "  npoints = " << npoints << std::endl;
  os << "  grid    = " << grid << std::endl;
  os << "  gdims   = " << gdims << std::endl;
  os << "  npair   = " << npair << std::endl;
  return true;
}

void VectorPairCorr::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = collectables.current();
  std::vector<RealType> tmp(npair*npoints);
  collectables.add(tmp.begin(), tmp.end());
}

void VectorPairCorr::registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  // add lattice information
  Matrix<double> axes(ndim);
  for (int l=0;l<ndim;l++)
    for (int m=0;m<ndim;m++)
      axes(l,m) = lattice.R(l,m);
  file.push(hdf_path{name_});
  file.write(axes, "axes");
  file.pop();
  std::vector<int> mesh(ndim);
  for (int l=0;l<ndim;l++) mesh[l] = grid[l];
  file.push(hdf_path{name_});
  file.write(mesh, "mesh");
  file.pop();
  // add data grid
  std::vector<int> ng(ndim+1);
  ng[0] = npair;
  for (int l=0;l<ndim;l++) ng[l+1] = grid[l];
  h5desc.emplace_back(hdf_path{name_});
  auto& oh = h5desc.back();
  oh.set_dimensions(ng, my_index_);
}

VectorPairCorr::Return_t VectorPairCorr::evaluate(ParticleSet& P)
{
  RealType wt = t_walker_->Weight;
  const auto& dii(P.getDistTableAA(d_aa_ID_));
  for (int iat = 1; iat < dii.centers(); ++iat)
  {
    const auto& drs = dii.getDisplRow(iat);
    const int ig = P.GroupID[iat];
    for (int j = 0; j < iat; ++j)
    { // steal from SpinDensity
      const int jg = P.GroupID[j];
      const int ipair = gen_pair_id(ig, jg, species.size());
      int offset = my_index_+ipair*npoints;
      // locate grid point
      const PosType u = lattice.toUnit(drs[j]);
      int point = 0;
      for (int l=0;l<ndim;l++)
        point += gdims[l] * ((int)(grid[l] * (u[l] - std::floor(u[l]))));
      // increament grid point
      if ((0 <= point) and (point < npoints)) P.Collectables[point+offset] += wt*norms(ig,jg);
    }
  }

  value_ = 0.0; // Value is no longer used in scalar.dat
  return value_;
}

int VectorPairCorr::gen_pair_id(const int ig, const int jg, const int ns) const
{ // steal from PairCorr
  if (jg < ig)
    return ns * (ns - 1) / 2 - (ns - jg) * (ns - jg - 1) / 2 + ig;
  else
    return ns * (ns - 1) / 2 - (ns - ig) * (ns - ig - 1) / 2 + jg;
}

std::unique_ptr<OperatorBase> VectorPairCorr::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<VectorPairCorr>(*this);
}

} // namespace qmcplusplus
