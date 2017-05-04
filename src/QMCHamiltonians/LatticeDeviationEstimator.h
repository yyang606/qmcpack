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
#ifndef QMCPLUSPLUS_LATTICEDEVIATION_H
#define QMCPLUSPLUS_LATTICEDEVIATION_H

#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "ParticleBase/ParticleAttribOps.h"
#include <Particle/DistanceTable.h>
#include <Particle/DistanceTableData.h>

namespace qmcplusplus
{

/** lattice deviation estimator
 *
 * Compute deviation of species="tgroup" in target particle set from species="sgroup" in source particle set. The motivation is to observe the deviation of protons from their crystal sites in an all electron-proton simulation of hydrogen i.e. two-component system
 */
class LatticeDeviationEstimator: public QMCHamiltonianBase
{
public:

  LatticeDeviationEstimator(ParticleSet& P, ParticleSet& sP, std::string tgroup, std::string sgroup);
  ~LatticeDeviationEstimator() { }

  bool put(xmlNodePtr cur);         // read input xml node, required
  bool get(std::ostream& os) const; // class description, required

  Return_t evaluate(ParticleSet& P); // main function that calculates the observable
  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  { // delegate responsity inline for speed
    return evaluate(P);
  }

  // allow multiple scalars to be registered in scalar.dat
  void addObservables(PropertySetType& plist, BufferType& collectables);
  void setObservables(PropertySetType& plist);

  // pure virtual functions require overrider
  void resetTargetParticleSet(ParticleSet& P);                            // required
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi); // required

private:
  SpeciesSet&       tspecies; // species table of target particle set
  SpeciesSet&       sspecies; // species table of source particle set
  std::string  tgroup,sgroup; // name of species to track
  DistanceTableData* d_table; // distance table between target and source particle sets
  ParticleSet&   tpset,spset; // save references to source and target particle sets
  int num_sites; // number of lattice sites (i.e. number of source particles)
  bool per_xyz;  // track deviation in each of x,y,z directions
  std::vector<RealType> xyz2; // temporary storage for deviation in each of x,y,z directions
  xmlNodePtr input_xml; // original xml
}; // LatticeDeviationEstimator

} // namespace qmcplusplus
#endif

/***************************************************************************
 * $RCSfile$   $Author: yyang $
 * $Revision: 7049 $   $Date: 2016-08-04 11:26:23 -0500 (Thur, 4 Aug 2017) $
 * $Id: LatticeDeviationEstimator.h 7049 2017-08-04 11:26:23 yyang $
 ***************************************************************************/

