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

#include "ComplexPolarization.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
ComplexPolarization::ComplexPolarization(ParticleSet& P) :
  tpset(P),
  lattice(P.getLattice()),
  ndim(lattice.ndim)
{
  // [v0_real, v0_imag]
  values.resize(2);
  // look along x by default
  axis = 0.0;
  axis[0] = 1.0;
  // maximum projcetion distance
  rmax = 0.0;
};

bool ComplexPolarization::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attrib;
  attrib.add(axis, "axis");
  attrib.put(cur);
  RealType norm=0.0;
  for (int l=0;l<axis.size();l++)
    norm += axis[l]*axis[l];
  axis /= std::sqrt(norm);
  auto axes = lattice.R;
  std::vector<RealType> lengths;
  for (int i=0;i<ndim;i++)
  {
    auto ai = axes.getRow(i);
    lengths.push_back(dot(ai, axis));
    for (int j=i+1;j<ndim;j++)
    {
      auto aj = axes.getRow(j);
      lengths.push_back(dot(ai+aj, axis));
    }
  }
  if (ndim == 3)
  {
    auto rvec = axes.getRow(0)+axes.getRow(1)+axes.getRow(2);
    lengths.push_back(dot(rvec, axis));
  }
  rmax = *std::max_element(lengths.begin(), lengths.end());
  get(app_log());
  return true;
}

bool ComplexPolarization::get(std::ostream& os) const
{ // class description
  os << "ComplexPolarization: " << name_ << std::endl;
  os << "  axis =" << axis << std::endl;
  os << "  rmax =" << rmax << std::endl;
  return true;
}

void ComplexPolarization::addObservables(PropertySetType& plist, BufferType& collectables)
{
  my_index_ = plist.size();
  for (int l=0;l<values.size();l++)
  { // make columns in scalar.dat
    std::ostringstream obsName;
    obsName << name_ << "_" << l;
    plist.add(obsName.str());
  }
}

void ComplexPolarization::setObservables(PropertySetType& plist)
{ // slots in plist must be allocated by addObservables() first
  copy(values.begin(), values.end(), plist.begin() + my_index_);
}

void ComplexPolarization::setParticlePropertyList(PropertySetType& plist, int offset)
{
  int index = my_index_ + offset;
  for (int i=0; i<values.size(); i++)
  {
    plist[index] = values[i];
    index++;
  }
}

ComplexPolarization::Return_t ComplexPolarization::evaluate(ParticleSet& P)
{
  RealType wgt = t_walker_->Weight;

  const int npart=P.getTotalNum();
  RealType cosz, sinz;

  // v = exp(i 2*pi * \sum_j f_j)
  RealType expo = 0.0;
  for (int j=0; j<npart; j++)
  {
    auto r = P.R[j];
    // fractional coordinate
    auto rp = dot(r, axis);
    auto fj = rp/rmax;
    expo += fj;
  }
  sincos(2*M_PI*expo, &cosz, &sinz);
  values[0] = cosz;
  values[1] = sinz;

  value_ = 0.0; // Value is no longer used in scalar.dat
  return value_;
}

std::unique_ptr<OperatorBase> ComplexPolarization::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<ComplexPolarization>(*this);
}

} // namespace qmcplusplus
