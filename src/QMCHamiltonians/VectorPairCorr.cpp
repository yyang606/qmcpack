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

namespace qmcplusplus
{
VectorPairCorr::VectorPairCorr(ParticleSet& P) :
  tpset(P),
  lattice(P.getLattice()),
  ndim(lattice.ndim),
  center(0.0),
  grid(-1)
{
  update_mode_.set(COLLECTABLE, 1);
};

bool VectorPairCorr::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attrib;
  attrib.put(cur);
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
      else if (name == "corner")
      {
        putContent(corner, element);
      }
    }
    element = element->next;
  }
  // number of grid points
  if (grid[0] <= 0)
  {
    std::ostringstream msg;
    msg << "grid is required in VectorPairCorr\n";
    throw std::runtime_error(msg.str());
  }
  npoints = 1;
  for (int l = 0; l < ndim; ++l)
    npoints *= grid[l];
  // corner
  corner = 0.0;
  for (int l = 0; l < ndim; ++l)
    corner[l] = center[l] - lattice.Center[l];
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
  os << "  corner  = " << corner << std::endl;
  return true;
}

void VectorPairCorr::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = collectables.current();
  std::vector<RealType> tmp(npoints);
  collectables.add(tmp.begin(), tmp.end());
}

void VectorPairCorr::registerCollectables(std::vector<ObservableHelper>& h5desc, hdf_archive& file) const
{
  std::vector<int> ng(1);
  ng[0] = npoints;

  hdf_path grp_name{name_};
  h5desc.emplace_back(grp_name);
  auto& oh = h5desc.back();
  oh.set_dimensions(ng, my_index_);
}

VectorPairCorr::Return_t VectorPairCorr::evaluate(ParticleSet& P)
{
  RealType wgt = t_walker_->Weight;
  value_ = 0.0; // Value is no longer used in scalar.dat
  return value_;
}

std::unique_ptr<OperatorBase> VectorPairCorr::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<VectorPairCorr>(*this);
}

} // namespace qmcplusplus
