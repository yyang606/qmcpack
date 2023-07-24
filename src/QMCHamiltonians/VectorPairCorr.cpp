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
  tpset(P),
  lattice(P.getLattice()),
  ndim(lattice.ndim),
  d_aa_ID_(P.addTable(P, DTModes::NEED_FULL_TABLE_ON_HOST_AFTER_DONEPBYP)),
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
  h5desc.emplace_back(hdf_path{name_});
  auto& oh = h5desc.back();
  oh.set_dimensions(ng, my_index_);
}

VectorPairCorr::Return_t VectorPairCorr::evaluate(ParticleSet& P)
{
  RealType wgt = t_walker_->Weight;
  const auto& dii(P.getDistTableAA(d_aa_ID_));
  int offset = my_index_;
  for (int iat = 1; iat < dii.centers(); ++iat)
  {
    const auto& drs = dii.getDisplRow(iat);
    for (int j = 0; j < iat; ++j)
    { // steal from SpinDensity
      const PosType u = lattice.toUnit(drs[j]);
      int point = offset;
      for (int l=0;l<ndim;l++)
        point += gdims[l] * ((int)(grid[l] * (u[l] - std::floor(u[l]))));
      if ((0 <= point-offset) and (point-offset < npoints)) P.Collectables[point] += wgt;
    }
  }

  value_ = 0.0; // Value is no longer used in scalar.dat
  return value_;
}

std::unique_ptr<OperatorBase> VectorPairCorr::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<VectorPairCorr>(*this);
}

} // namespace qmcplusplus
