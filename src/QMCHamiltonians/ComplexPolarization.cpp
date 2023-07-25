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

namespace qmcplusplus
{
ComplexPolarization::ComplexPolarization(ParticleSet& P) :
  tpset(P),
  lattice(P.getLattice()),
  ndim(lattice.ndim)
{
  // [v0_real, v0_imag, v1_real, v1_imag, ...]
  values.resize(2*ndim);
};

bool ComplexPolarization::put(xmlNodePtr cur)
{
  return true;
}

bool ComplexPolarization::get(std::ostream& os) const
{ // class description
  os << "ComplexPolarization: " << name_;
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
  auto lattice = P.getLattice();
  RealType cosz, sinz;

  std::vector<RealType> expos(ndim);
  std::fill(expos.begin(), expos.end(), 0.0);

  // v_l = exp(i 2*pi * \sum_j f_{jl})
  for (int i=0; i<npart; i++)
  {
    auto r = P.R[i];
    // fractional coordinate
    auto f = lattice.toUnit(r);
    for (int l=0; l<ndim; l++) expos[l] += f[l];
  }
  for (int l=0; l<ndim; l++)
  {
    sincos(2*M_PI*expos[l], &cosz, &sinz);
    values[2*l] = cosz;
    values[2*l+1] = sinz;
  }

  value_ = 0.0; // Value is no longer used in scalar.dat
  return value_;
}

std::unique_ptr<OperatorBase> ComplexPolarization::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  return std::make_unique<ComplexPolarization>(*this);
}

} // namespace qmcplusplus
