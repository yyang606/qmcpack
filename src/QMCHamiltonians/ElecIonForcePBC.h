//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Yubo Yang, paul.young.0414@gmail.com, University of Illinois at Urbana-Champaign
//  mostly copied from ForceChiesaPBCAA
//
// File created by: Yubo Yang, paul.young.0414@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ELEC_ION_FORCE_PBC_H
#define QMCPLUSPLUS_ELEC_ION_FORCE_PBC_H

#include <QMCHamiltonians/QMCHamiltonianBase.h>
#include <LongRange/LRCoulombSingleton.h>

namespace qmcplusplus
{

class ElecIonForcePBC: public QMCHamiltonianBase
{
public:
  typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
  typedef ParticleSet::ParticlePos_t        ForceStorageType;

  ElecIonForcePBC(ParticleSet& P, ParticleSet& sP);
  ~ElecIonForcePBC(){};

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

  // make room in hdf5 observable registry
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const;

  // pure virtual functions require overrider
  void resetTargetParticleSet(ParticleSet& P);                            // required
  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi); // required

private:
  // I/O info
  bool hdf5_out;        // use .h5 file for data (follow SkEstimator)
  int  h5_index;        // track the starting memory location in P.Collectables
  xmlNodePtr input_xml; // pointer to xml <estimator> entry

  // particle info
  SpeciesSet&       tspecies; // species table of target particle set
  SpeciesSet&       sspecies; // species table of source particle set
  ParticleSet&   tpset,spset; // save references to source and target particle sets
  int nelec,nion,nespec,nispec;
  std::vector<RealType> elec_charge_vec,ion_charge_vec;

  int tid; // index of electron-ion distance table

  // force estimator parameters
  ForceStorageType forces; // array to store forces in
  RealType Rcut;         // radial distance within which estimator is used
  int m_exp;             // exponent in polynomial fit
  int N_basis;           // size of polynomial basis
  Matrix<RealType> Sinv; // terms in fitting polynomial
  Vector<RealType> h;    // terms in fitting polynomial
  Vector<RealType> c;    // polynomial coefficients

  // local methods
  RealType g_filter(RealType r);
  void initSRFit(RealType Rcut,int m_exp,int N_basis,
    Matrix<RealType> &Sinv, Vector<RealType> &h, Vector<RealType> &c);

  // short-range forces
  void evaluateSR(ParticleSet& P,ForceStorageType& forces);

  // long-range Coulomb stuff
  LRHandlerType* AB;
  void evaluateLR(ParticleSet& P,ForceStorageType& forces);

}; // ElecIonForcePBC

} // namespace qmcplusplus
#endif
