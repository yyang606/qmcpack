//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCHamiltonians/SpeciesKineticEnergy.h"
#include "QMCHamiltonians/BareKineticEnergy.h" // laplaician() defined here
#include <OhmmsData/AttributeSet.h>

namespace qmcplusplus
{

SpeciesKineticEnergy::SpeciesKineticEnergy(ParticleSet& P) : tpset(P), hdf5_out(false)
{
  SpeciesSet& tspecies(P.getSpeciesSet());
  int massind = tspecies.getAttribute("mass");

  num_species = tspecies.size();
  app_log() << "SpeciesKineticEnergy is tracking " << num_species << " species." << std::endl;

  species_kinetic.resize(num_species);
  vec_minus_over_2m.resize(num_species);
  species_names.resize(num_species);
  for (int ispec=0; ispec<num_species; ispec++)
  {
    species_names[ispec] = tspecies.speciesName[ispec];
    vec_minus_over_2m[ispec] = -1./(2.*tspecies(massind,ispec)); 
    app_log() << "  " << species_names[ispec] << " mass = " << vec_minus_over_2m[ispec] << std::endl;
  }

};

bool SpeciesKineticEnergy::put(xmlNodePtr cur)
{ 
  // read hdf5="yes" or "no"
  OhmmsAttributeSet attrib;
  std::string hdf5_flag="no";
  std::string scalar_flag="yes";
  attrib.add(hdf5_flag,"hdf5");
  attrib.put(cur);

  if (hdf5_flag=="yes")
    hdf5_out=true;
  else if (hdf5_flag=="no")
    hdf5_out=false;
  else
    APP_ABORT("SpeciesKineticEnergy::put() - Please choose \"yes\" or \"no\" for hdf5 flag");
  // end if

  if (hdf5_out)
  { // add this estimator to stat.h5
    UpdateMode.set(COLLECTABLE,1);
  }
  return true;
}

bool SpeciesKineticEnergy::get(std::ostream& os) const
{ // class description
  os << "SpeciesKineticEnergy: " << myName;
  return true;
}

void SpeciesKineticEnergy::addObservables(PropertySetType& plist, BufferType& collectables)
{
  myIndex = plist.size();
  for (int ispec=0; ispec<num_species; ispec++)
  { // make columns named "$myName_u", "$myName_d" etc.
    plist.add(myName + "_" + species_names[ispec]);
  }
  if (hdf5_out)
  { // make room in h5 file and save its index
    h5_index = collectables.size();
    std::vector<RealType> tmp(num_species);
    collectables.add(tmp.begin(),tmp.end());
  }
}

void SpeciesKineticEnergy::registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const
{
  if (hdf5_out)
  {
    std::vector<int> ndim(1,num_species);
    observable_helper* h5o=new observable_helper(myName);
    h5o->set_dimensions(ndim,h5_index);
    h5o->open(gid);
    h5desc.push_back(h5o);
  }
}

void SpeciesKineticEnergy::setObservables(PropertySetType& plist)
{ // slots in plist must be allocated by addObservables() first
  copy(species_kinetic.begin(),species_kinetic.end(),plist.begin()+myIndex);
}

SpeciesKineticEnergy::Return_t SpeciesKineticEnergy::evaluate(ParticleSet& P)
{
  std::fill(species_kinetic.begin(),species_kinetic.end(),0.0);
  RealType wgt = tWalker->Weight;

  for (int iat=0; iat<P.getTotalNum(); iat++)
  {
    int ispec  = P.GroupID(iat);
    RealType my_kinetic = vec_minus_over_2m[ispec]*laplacian(P.G[iat],P.L[iat]);
    if (hdf5_out)
    {
      P.Collectables[h5_index + ispec] += my_kinetic*wgt;
    }
    species_kinetic[ispec] += my_kinetic;
  }
  
  Value = 0.0; // Value is no longer used
  return Value;
}

QMCHamiltonianBase* SpeciesKineticEnergy::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{ 
  return new SpeciesKineticEnergy(*this);
}

} // namespace qmcplusplus

/***************************************************************************
 * $RCSfile$   $Author: yyang $
 * $Revision: 7049 $   $Date: 2016-08-03 20:47:23 -0500 (Wed, 3 Aug 2017) $
 * $Id: SpeciesKineticEnergy.h 7049 2017-08-03 20:47:23 yyang $
 ***************************************************************************/

