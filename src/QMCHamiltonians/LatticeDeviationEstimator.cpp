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
#include "QMCHamiltonians/LatticeDeviationEstimator.h"
#include <OhmmsData/AttributeSet.h>

namespace qmcplusplus
{

LatticeDeviationEstimator::LatticeDeviationEstimator(ParticleSet& P,ParticleSet& sP, std::string tgroup_in, std::string sgroup_in): 
  tspecies(P.getSpeciesSet()), sspecies(sP.getSpeciesSet()), spset(sP), tpset(P)
{ 
  // get the distance table from quantum particle set
  int tid   = P.getTable(sP);
  d_table   = P.DistTables[tid];
  tgroup    = tgroup_in;
  sgroup    = sgroup_in;
  per_xyz   = false;
  
  // calculate number of source particles to use as lattice sites
  int src_species_id = sspecies.findSpecies(sgroup);
  num_sites = spset.last(src_species_id) - spset.first(src_species_id);
  int tar_species_id = tspecies.findSpecies(tgroup);
  int num_tars = tpset.last(tar_species_id) - tpset.first(tar_species_id);
  if (num_tars!=num_sites)
  {
    app_log() << "number of target particles = " << num_tars << std::endl;
    app_log() << "number of source particles = " << num_sites << std::endl;
    APP_ABORT("nsource != ntargets in LatticeDeviationEstimator");
  }
}

bool LatticeDeviationEstimator::put(xmlNodePtr cur)
{ 
  input_xml = cur;
  std::string xyz_flag="no";
  OhmmsAttributeSet attrib;
  attrib.add(xyz_flag,"per_xyz");
  attrib.put(cur);

  if (xyz_flag=="yes")
  {
    per_xyz = true;
    xyz2.resize(OHMMS_DIM);
  } else if (xyz_flag=="no") {
    per_xyz = false;
  } else {
    APP_ABORT("LatticeDeviationEstimator::put() - Please choose \"yes\" or \"no\" for per_xyz flag");
  } // end if xyz_flag
  
  return true;
}

bool LatticeDeviationEstimator::get(std::ostream& os) const
{ // class description
  os << "LatticeDeviationEstimator: " << myName << "lattice = " << spset.getName();
  return true;
}

LatticeDeviationEstimator::Return_t LatticeDeviationEstimator::evaluate(ParticleSet& P)
{ // calculate <r^2> averaged over lattice sites
  Value = 0.0;
  std::fill(xyz2.begin(),xyz2.end(),0.0);

  RealType wgt = tWalker->Weight;

  // temp variables
  RealType r2;
  PosType  dr;
  
  int nsite(0);    // site index
  int cur_jat(-1); // target particle index
  for (int iat=0;iat<spset.getTotalNum();iat++)
  { // for each desired source particle
    if (sspecies.speciesName[spset.GroupID(iat)] == sgroup)
    { // find desired species
      for (int jat=cur_jat+1;jat<tpset.getTotalNum();jat++)
      { // find corresponding (!!!! assume next) target particle
        if (tspecies.speciesName[tpset.GroupID(jat)] == tgroup)
        {
          // distance between particle iat in source pset, and jat in target pset
          int nn = d_table->loc(iat,jat); // location where distance is stored 
          r2 = std::pow( d_table->r( nn ) ,2);
          Value += r2;

          if (per_xyz)
          {
            dr = d_table->dr(nn);
            for (int idir=0; idir<OHMMS_DIM; idir++)
            {
              RealType dir2 = dr[idir]*dr[idir]; 
              xyz2[idir] += dir2;
            }
          }

          cur_jat = jat;
          break;
        }
      }
      nsite += 1; // count the number of sites, for checking only
    } // found desired species (source particle)
  }

  if (nsite!=num_sites) 
  {
    app_log() << "num_sites = " << num_sites << " nsites = " << nsite << std::endl;
    APP_ABORT("Mismatch in LatticeDeivationEstimator.");
  }

  // average per site
  Value /= num_sites;
  if (per_xyz)
  {
    std::transform(xyz2.begin(),xyz2.end(),xyz2.begin(),
        bind2nd(std::multiplies<RealType>(),1./num_sites) );
  }

  return Value;
}

void LatticeDeviationEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if (per_xyz)
  {
    myIndex = plist.size();
    for (int idir=0; idir<OHMMS_DIM; idir++)
    {
      std::stringstream ss;
      ss << idir;
      plist.add(myName + "_dir" + ss.str() );
    }
  } else {
    myIndex = plist.add(myName); // same as QMCHamiltonianBase::addObservables
  }
}

void LatticeDeviationEstimator::setObservables(PropertySetType& plist)
{
  if (per_xyz)
  {
    copy(xyz2.begin(),xyz2.end(),plist.begin()+myIndex);
  } else {
    plist[myIndex] = Value; // default behavior
  }
}

void LatticeDeviationEstimator::resetTargetParticleSet(ParticleSet& P)
{
}

QMCHamiltonianBase* LatticeDeviationEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{ 
  LatticeDeviationEstimator* myclone = new LatticeDeviationEstimator(*this); 
  //myclone->tpset = &qp;

  int tid = qp.getTable(spset);
  myclone->d_table = qp.DistTables[tid];
  
  return myclone;
}

} // namespace qmcplusplus

/***************************************************************************
 * $RCSfile$   $Author: yyang $
 * $Revision: 7049 $   $Date: 2016-08-04 11:26:23 -0500 (Thur, 4 Aug 2017) $
 * $Id: LatticeDeviationEstimator.h 7049 2017-08-04 11:26:23 yyang $
 ***************************************************************************/

