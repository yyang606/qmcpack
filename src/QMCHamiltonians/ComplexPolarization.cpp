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
  // look along b_1 by default
  gvec = 0;
  gvec[0] = 1;
};

bool ComplexPolarization::put(xmlNodePtr cur)
{
  OhmmsAttributeSet attrib;
  attrib.add(gvec, "gvec");
  attrib.put(cur);
  get(app_log());
  return true;
}

bool ComplexPolarization::get(std::ostream& os) const
{ // class description
  os << "ComplexPolarization: " << name_ << std::endl;
  os << "  gvec =" << gvec << std::endl;
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
    auto rj = P.R[j];
    // fractional coordinate
    auto sj = lattice.toUnit(rj);
    auto grj = 2.0*M_PI*dot(gvec, sj);
    expo += grj;
  }
  sincos(expo, &cosz, &sinz);
  values[0] = wgt*cosz;
  values[1] = wgt*sinz;

  value_ = 0.0; // Value is no longer used in scalar.dat
  return value_;
}

std::unique_ptr<OperatorBase> ComplexPolarization::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<ComplexPolarization>(*this);
}

} // namespace qmcplusplus
