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

#include <QMCHamiltonians/ElecIonForcePBC.h>
#include <Particle/DistanceTableData.h>
#include <OhmmsData/AttributeSet.h>
#include <OhmmsData/ParameterSet.h>
#include <Numerics/DeterminantOperators.h>
#include <Numerics/MatrixOperators.h>

namespace qmcplusplus
{

ElecIonForcePBC::ElecIonForcePBC(ParticleSet& P,ParticleSet& sP):
  spset(sP), tpset(P),
  tspecies(P.getSpeciesSet()), sspecies(sP.getSpeciesSet()), 
  hdf5_out(false),
  Rcut(0.4), N_basis(2), m_exp(2)
{
  // get electron-ion distance table, create if not already initialized
  tid   = P.addTable(sP,DT_AOS); // getTable(sP) does not work with threads

  // save particle info.
  nelec =  P.getTotalNum(); // target
  nion  = sP.getTotalNum(); // source
  nespec= tspecies.size();  // number of target species (probably 2 for "u","d")
  nispec= sspecies.size();  // number of source species
  // get particle charges
  elec_charge_vec.resize(nelec);
  int elec_charge_idx = tspecies.getAttribute("charge");
  for (int ielec=0;ielec<nelec;ielec++)
  {
    int is = P.GroupID(ielec); // species ID
    elec_charge_vec[ielec] = tspecies(elec_charge_idx,is);
  }
  ion_charge_vec.resize(nion);
  int ion_charge_idx = sspecies.getAttribute("charge");
  for (int iat=0;iat<nion;iat++)
  {
    int is = sP.GroupID(iat);
    ion_charge_vec[iat] = sspecies(ion_charge_idx,is);
  }
  
  // store as many force vectors as there are source (classical) particles
  forces.resize(nion);
  forces = 0.0;

  // initialize long-range (LR) breakup
  AB = LRCoulombSingleton::getDerivHandler(P);

  // initSRFit must be called (in put()) to enable short-range (SR) evaluation
}

bool ElecIonForcePBC::put(xmlNodePtr cur)
{
  // save xml reference for myclone
  input_xml = cur;

  // parse xml attributes
  std::string hdf5_in="no";
  OhmmsAttributeSet attrib;
  attrib.add(hdf5_in,"hdf5");
  attrib.put(cur);

  if (hdf5_in=="yes")
  {
    hdf5_out = true;
  } else if (hdf5_in=="no") 
  {
    hdf5_out = false;
  } else 
  {
    APP_ABORT("ElecIonForcePBC::put() - Please choose \"yes\" or \"no\" for hdf5");
  } // end if hdf5_in

  if (hdf5_out) UpdateMode.set(COLLECTABLE,1);

  // parse xml <parameter>
  ParameterSet force_params;
  force_params.add(Rcut,"rcut","real");
  force_params.add(N_basis,"nbasis","int");
  force_params.add(m_exp,"weight_exp","int");
  force_params.put(cur);

  this->get(app_log()); // output a string rep. of this estimator
  // useful for checking that input values are accepted

  // allocate short-range matrix
  Sinv.resize(N_basis, N_basis);
  h.resize(N_basis);
  c.resize(N_basis);

  initSRFit(Rcut,m_exp,N_basis // inputs
    ,Sinv,h,c // outputs
  );
}

bool ElecIonForcePBC::get(std::ostream& os) const
{
  os << "ElecIonForcePBC Parameters:" << std::endl;
  os << "  Rcut    = " << Rcut        << std::endl;
  os << "  N_basis = " << N_basis     << std::endl;
  os << "  m_exp   = " << m_exp       << std::endl;
  os << "  hdf5    = " << hdf5_out    << std::endl;
}

ElecIonForcePBC::Return_t ElecIonForcePBC::evaluate(ParticleSet& P)
{
  RealType ftot = 0.0;
  forces = 0.0;
  evaluateLR(P,forces); // add short-range contributions to forces
  evaluateSR(P,forces); // add long-range contributions to forces

  // get total force, fill P.Collectables if needed
  for (int iat=0;iat<nion;iat++)
  {
    for (int idim=0;idim<OHMMS_DIM;idim++)
    {
      RealType myf = forces[iat][idim];
      ftot += myf;
      if (hdf5_out)
      {
        RealType wgt = tWalker->Weight;
        P.Collectables[h5_index+iat*OHMMS_DIM+idim] += wgt*myf;
      }
    }
  }

  Value = ftot;
  return Value;
}

void ElecIonForcePBC::addObservables(PropertySetType& plist, BufferType& collectables)
{
  // get myIndex for scalar.dat
  myIndex = plist.size();
  if (!hdf5_out)
  { // use scalar.dat only
    for (int iat=0;iat<nion;iat++)
    {
      for (int idim=0;idim<OHMMS_DIM;idim++)
      {
        std::stringstream obsName;
        obsName << myName << "_" << iat << "_" << idim;
        plist.add(obsName.str());
      }
    }
  } else
  { // use stat.h5
    h5_index = collectables.size();
    std::vector<RealType> tmp;
    tmp.resize(nion*OHMMS_DIM);
    collectables.add(tmp.begin(),tmp.end());
    // store total forces in scalar.dat
    std::stringstream obsName;
    obsName << myName;
    plist.add(obsName.str());
  }
}

void ElecIonForcePBC::setObservables(PropertySetType& plist)
{
  if (!hdf5_out)
  { // store particle-resolved force in stat.h5
    int offset = myIndex;
    for (int iat=0;iat<nion;iat++)
    {
      for (int idim=0;idim<OHMMS_DIM;idim++)
      {
        plist[offset] = forces[iat][idim];
        offset++;
      }
    }
  } else
  { // store total force in scalar.dat
    plist[myIndex] = Value;
  }
}

void ElecIonForcePBC::registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const
{
  if (hdf5_out)
  {
    std::vector<int> ndim(1,nion*OHMMS_DIM);

    // open hdf5 entry and resize
    observable_helper* h5o=new observable_helper(myName);
    h5o->set_dimensions(ndim,h5_index);
    h5o->open(gid);
  
    // add to h5 file
    h5desc.push_back(h5o);
  }
}

void ElecIonForcePBC::resetTargetParticleSet(ParticleSet& P)
{
  int new_tid = P.addTable(spset,DT_AOS);
  if (new_tid != tid)
  {
    APP_ABORT("ElecIonForcePBC found inconsistent distance table index");
  }
  AB->resetTargetParticleSet(P);
}

QMCHamiltonianBase* ElecIonForcePBC::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  ElecIonForcePBC* myclone = new ElecIonForcePBC(qp,spset);
  myclone->put(input_xml);
  return myclone;
}

void ElecIonForcePBC::initSRFit(RealType Rcut,int m_exp,int N_basis,
  Matrix<RealType> &Sinv, Vector<RealType> &h, Vector<RealType> &c)
{
  for(int k=0; k<N_basis; k++)
  {
    h[k] = std::pow(Rcut, (k+2))/static_cast<RealType>(k+2);
    for(int j=0; j<N_basis; j++)
    {
      Sinv(k,j) = std::pow(Rcut, (m_exp+k+j+3))/static_cast<RealType>(m_exp+k+j+3);
    }
  }
  // in Numerics/DeterminantOperators.h
  invert_matrix(Sinv, false);
  // in Numerics/MatrixOperators.h
  MatrixOperators::product(Sinv, h.data(), c.data());
}


ElecIonForcePBC::RealType ElecIonForcePBC::g_filter(RealType r)
{
  RealType g_q=0.0;
  for (int q=0; q<N_basis; q++)
  {
   g_q += c[q]*std::pow(r,m_exp+q+1);
  }
  return g_q;
}

void ElecIonForcePBC::evaluateSR(ParticleSet& P,ForceStorageType& forces)
{
  const DistanceTableData &d_ab(*P.DistTables[tid]);
  for (int iat=0;iat<nion;iat++)
  {
    for (int nn=d_ab.M[iat],jat=0;nn<d_ab.M[iat+1];++nn,++jat)
    {
      RealType rij = d_ab.r(nn), inv_rij = d_ab.rinv(nn);
      RealType g_f = (rij>=Rcut) ? 1.0 : g_filter(rij);
      RealType V = -AB->srDf(rij,inv_rij);
      PosType drhat = inv_rij * d_ab.dr(nn); // direction

      forces[iat] += -g_f*ion_charge_vec[iat]*elec_charge_vec[jat]*V*drhat;
    }
  }
}
void ElecIonForcePBC::evaluateLR(ParticleSet& P,ForceStorageType& forces)
{
  std::vector<TinyVector<RealType,DIM> > grad(nion);
  for (int j=0; j<nespec; j++)
  {
    for (int iat=0; iat<nion; iat++)
      grad[iat] = TinyVector<RealType,DIM>(0.0, 0.0, 0.0);
    AB->evaluateGrad(spset, P, j, ion_charge_vec, grad);
    for (int iat=0; iat<grad.size(); iat++){
      forces[iat] += elec_charge_vec[j]*grad[iat];
    }
  } // end for electron species
}

} // namespace qmcplusplus
