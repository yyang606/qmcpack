//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Yubo Yang, paul.young.0414@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Yubo Yang, paul.young.0414@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "LatticeDeviationEstimator.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
LatticeDeviationEstimator::LatticeDeviationEstimator(ParticleSet& P,
                                                     ParticleSet& sP,
                                                     const std::string& tgroup_in,
                                                     const std::string& sgroup_in)
    : tspecies(P.getSpeciesSet()),
      sspecies(sP.getSpeciesSet()),
      tpset(P),
      spset(sP),
      tgroup(tgroup_in),
      sgroup(sgroup_in),
      hdf5_out(false),
      per_xyz(false),
      myTableID_(P.addTable(sP)),
      iat_first(-1),
      iat_last(-1),
      iel_first(-1),
      iel_last(-1)
{
  // calculate number of source particles to use as lattice sites
  int src_species_id = sspecies.findSpecies(sgroup);
  if (src_species_id == sspecies.size())
    throw std::runtime_error("source group not found");
  iat_first = spset.first(src_species_id);
  iat_last = spset.last(src_species_id);
  num_sites = iat_last-iat_first;
  int tar_species_id = tspecies.findSpecies(tgroup);
  if (tar_species_id == tspecies.size())
    throw std::runtime_error("target group not found");
  iel_first = tpset.first(tar_species_id);
  iel_last = tpset.last(tar_species_id);
  int num_tars = iel_last-iel_first;
  if (num_tars != num_sites)
  {
    app_log() << "number of target particles = " << num_tars << std::endl;
    app_log() << "number of source particles = " << num_sites << std::endl;
    APP_ABORT("nsource != ntargets in LatticeDeviationEstimator");
  }
  rij.resize(num_tars);
  for (int i=0;i<num_tars;i++)
  {
    rij[i].resize(num_sites);
  }
  icols.resize(num_sites);
}

bool LatticeDeviationEstimator::put(xmlNodePtr cur)
{
  input_xml             = cur;
  std::string hdf5_flag = "no";
  std::string xyz_flag  = "no";
  OhmmsAttributeSet attrib;
  attrib.add(hdf5_flag, "hdf5");
  attrib.add(xyz_flag, "per_xyz");
  attrib.put(cur);

  if (hdf5_flag == "yes")
  {
    hdf5_out = true;
  }
  else if (hdf5_flag == "no")
  {
    hdf5_out = false;
  }
  else
  {
    APP_ABORT("LatticeDeviationEstimator::put() - Please choose \"yes\" or \"no\" for hdf5 flag");
  } // end if hdf5_flag
  if (hdf5_out)
  {
    // change the default behavior of addValue() called by addObservables()
    // YY: does this still matter if I have overwritten addObservables()?
    UpdateMode.set(COLLECTABLE, 1);
  }

  if (xyz_flag == "yes")
  {
    per_xyz = true;
    xyz2.resize(OHMMS_DIM);
  }
  else if (xyz_flag == "no")
  {
    per_xyz = false;
  }
  else
  {
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
  std::fill(xyz2.begin(), xyz2.end(), 0.0);

  RealType wgt        = tWalker->Weight;
  const auto& d_table = P.getDistTable(myTableID_);

  // temp variables
  RealType r, r2;
  PosType dr;
  // extract distance table
  for (int iat=iat_first,i=0; iat<iat_last; iat++,i++)
  {
      for (int jat=iel_first,j=0; jat<iel_last; jat++,j++)
      {
          r = d_table.getDistRow(jat)[iat];
          rij[i][j] = r;
      }
  }
  // linear sum assignment
  for (int i=0;i<icols.size();i++)
  {
    icols[i] = i; // !!!! place holder
  }
  // fill output
  for (int iat = iat_first; iat < iat_last; iat++)
  {
    r = rij[iat][icols[iat]];
    r2 = r * r;
    Value += r2;

    if (hdf5_out & !per_xyz)
    { // store deviration for each lattice site if h5 file is available
      P.Collectables[h5_index + num_sites] = wgt * r2;
    }

    if (per_xyz)
    {
      //dr = d_table.getDisplRow(jat)[iat]; // !!!! need fix
      for (int idir = 0; idir < OHMMS_DIM; idir++)
      {
        RealType dir2 = dr[idir] * dr[idir];
        xyz2[idir] += dir2;
        if (hdf5_out)
        {
          P.Collectables[h5_index + num_sites * OHMMS_DIM + idir] = wgt * dir2;
        }
      }
    }
  }

  // average per site
  Value /= num_sites;
  if (per_xyz)
  {
    for (int idir = 0; idir < OHMMS_DIM; idir++)
    {
      xyz2[idir] /= num_sites;
    }
  }

  return Value;
}

void LatticeDeviationEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  // get myIndex for scalar.dat
  if (per_xyz)
  {
    myIndex = plist.size();
    for (int idir = 0; idir < OHMMS_DIM; idir++)
    {
      std::stringstream ss;
      ss << idir;
      plist.add(myName + "_dir" + ss.str());
    }
  }
  else
  {
    myIndex = plist.add(myName); // same as OperatorBase::addObservables
  }

  // get h5_index for stat.h5
  if (hdf5_out)
  {
    h5_index = collectables.size();
    std::vector<RealType> tmp;
    if (per_xyz)
    {
      tmp.resize(num_sites * OHMMS_DIM);
    }
    else
    {
      tmp.resize(num_sites);
    }
    collectables.add(tmp.begin(), tmp.end());
  }
}

void LatticeDeviationEstimator::setObservables(PropertySetType& plist)
{
  if (per_xyz)
  {
    copy(xyz2.begin(), xyz2.end(), plist.begin() + myIndex);
  }
  else
  {
    plist[myIndex] = Value; // default behavior
  }
}

void LatticeDeviationEstimator::resetTargetParticleSet(ParticleSet& P) {}

OperatorBase* LatticeDeviationEstimator::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  // default constructor does not work with threads
  //LatticeDeviationEstimator* myclone = new LatticeDeviationEstimator(*this);
  LatticeDeviationEstimator* myclone = new LatticeDeviationEstimator(qp, spset, tgroup, sgroup);
  myclone->put(input_xml);

  return myclone;
}

void LatticeDeviationEstimator::registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const
{
  if (hdf5_out)
  {
    // one scalar per lattice site (i.e. source particle)
    std::vector<int> ndim(1, num_sites);
    if (per_xyz)
    {
      // one scalar per lattice site per dimension
      ndim[0] = num_sites * OHMMS_DIM;
    }

    // open hdf5 entry and resize
    observable_helper* h5o = new observable_helper(myName);
    h5o->set_dimensions(ndim, h5_index);
    h5o->open(gid);

    // add to h5 file
    h5desc.push_back(h5o);
  }
}


} // namespace qmcplusplus
