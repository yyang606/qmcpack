//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_DIFFERENTIAL_TWOBODYJASTROW_H
#define QMCPLUSPLUS_DIFFERENTIAL_TWOBODYJASTROW_H
#include "Configuration.h"
#include "QMCWaveFunctions/DiffOrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Utilities/IteratorUtility.h"

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 *  @brief Specialization for two-body Jastrow function using multiple functors
 */
template<class FT>
class DiffTwoBodyJastrowOrbital: public DiffOrbitalBase
{
  ///number of variables this object handles
  int NumVars;
  ///number of target particles
  int NumPtcls;
  ///number of groups, e.g., for the up/down electrons
  int NumGroups;
  ///variables handled by this orbital
  opt_variables_type myVars;
  ///container for the Jastrow functions  for all the pairs
  std::vector<FT*> F;
  ///offset for the optimizable variables
  std::vector<std::pair<int,int> > OffSet;
  Vector<RealType> dLogPsi;
  std::vector<GradVectorType*> gradLogPsi;
  std::vector<ValueVectorType*> lapLogPsi;
  std::map<std::string,FT*> J2Unique;
  std::vector<RealType> mass_vec;
  ParticleSet* PtclRef;

public:

  ///constructor
  DiffTwoBodyJastrowOrbital(ParticleSet& p):NumVars(0)
  {
    PtclRef = &p;
    NumPtcls=p.getTotalNum();
    NumGroups=p.groups();
    F.resize(NumGroups*NumGroups,0);
    
    // save particle masses for kinetic energy (needed in evaluateDerivatives)
    mass_vec.resize(NumPtcls);
    SpeciesSet& tspecies(p.getSpeciesSet());
    int massind = tspecies.getAttribute("mass");
    int num_species = tspecies.size();
    for (int iptcl=0;iptcl<NumPtcls;iptcl++)
    {
      int ispec = p.GroupID(iptcl);
      mass_vec[iptcl] = tspecies(massind,ispec);
    }

  }

  ~DiffTwoBodyJastrowOrbital()
  {
    delete_iter(gradLogPsi.begin(),gradLogPsi.end());
    delete_iter(lapLogPsi.begin(),lapLogPsi.end());
  }

  void addFunc(int ia, int ib, FT* rf)
  { // add radial function rf to J2Unique
    SpeciesSet& species(PtclRef->getSpeciesSet());
    std::string pair_name = species.speciesName[ia] + species.speciesName[ib];
    J2Unique[pair_name]=rf;

    // also assign rf to correlate species pairs (ia,ib), (ib,ia)
    F[ia*NumGroups+ib]=rf;
    F[ib*NumGroups+ia]=rf; // enforce exchange symmetry
  }

  void linkFunc(int ia, int ib, FT* rf, bool exchange=true)
  { // assign FunctorType rf to species pair (ia,ib)
    F[ia*NumGroups+ib] = rf;
    if (exchange) F[ib*NumGroups+ia] = rf;
  }

  FT* getFunc(int ia, int ib)
  { // locate pair function by species indices (ia,ib)
    return F[ia*NumGroups+ib];
  }

  FT* findFunc(std::string pair_name)
  { // locate pair function by species pair name such as 'uu'
    FT* rf;
    bool found_coeff = false;
    typename std::map<std::string,FT*>::const_iterator
      it(J2Unique.begin()), it_end(J2Unique.end());
    for (;it!=it_end;it++)
    { // it->first: 'uu','ud'
      if (it->first==pair_name)
      {
        rf = it->second;
        found_coeff = true;
      }
    }
    if (!found_coeff) APP_ABORT(pair_name+" not found");
    return rf;
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
    typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
    while (it != it_end)
    {
      (*it++).second->resetParameters(active);
    }
  }

  ///reset the distance table
  void resetTargetParticleSet(ParticleSet& P)
  {
    PtclRef = &P;
  }

  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.clear();
    typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
    while (it != it_end)
    {
      (*it).second->myVars.getIndex(active);
      myVars.insertFrom((*it).second->myVars);
      ++it;
    }
    myVars.getIndex(active);
    NumVars=myVars.size();

    //myVars.print(std::cout);

    if (NumVars && dLogPsi.size()==0)
    {
      dLogPsi.resize(NumVars);
      gradLogPsi.resize(NumVars,0);
      lapLogPsi.resize(NumVars,0);
      for (int i=0; i<NumVars; ++i)
      {
        gradLogPsi[i]=new GradVectorType(NumPtcls);
        lapLogPsi[i]=new ValueVectorType(NumPtcls);
      }
      OffSet.resize(F.size());
      int varoffset=myVars.Index[0];
      for (int i=0; i<F.size(); ++i)
      {
        if(F[i] && F[i]->myVars.Index.size())
        {
          OffSet[i].first=F[i]->myVars.Index.front()-varoffset;
          OffSet[i].second=F[i]->myVars.Index.size()+OffSet[i].first;
        }
        else
        {
          OffSet[i].first=OffSet[i].second=-1;
        }
      }
    }
  }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi)
  {
    if(myVars.size()==0) return;
    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(),false);
    for (int k=0; k<myVars.size(); ++k)
    {
      int kk=myVars.where(k);
      if (kk<0)
        continue;
      if (active.recompute(kk))
        recalculate=true;
      rcsingles[k]=true;
    }
    if (recalculate)
    {
      ///precomputed recalculation switch
      std::vector<bool> RecalcSwitch(F.size(),false);
      for (int i=0; i<F.size(); ++i)
      {
        if(OffSet[i].first<0)
        {
          // nothing to optimize
          RecalcSwitch[i]=false;
        }
        else
        {
          bool recalcFunc(false);
          for (int rcs=OffSet[i].first; rcs<OffSet[i].second; rcs++)
            if (rcsingles[rcs]==true) recalcFunc=true;
          RecalcSwitch[i]=recalcFunc;
        }
      }
      dLogPsi=0.0;
      for (int p=0; p<NumVars; ++p)
        (*gradLogPsi[p])=0.0;
      for (int p=0; p<NumVars; ++p)
        (*lapLogPsi[p])=0.0;
      std::vector<TinyVector<RealType,3> > derivs(NumVars);
      const DistanceTableData* d_table=P.DistTables[0];
      if(d_table->DTType == DT_SOA)
      {
        constexpr RealType cone(1);
        constexpr RealType lapfac(OHMMS_DIM-cone);
        const size_t n=d_table->size(SourceIndex);
        const size_t ng=P.groups();
        for(size_t i=1; i<n; ++i)
        {
          const size_t ig=P.GroupID[i]*ng;
          const RealType* dist=d_table->Distances[i];
          const auto& displ=d_table->Displacements[i];
          for(size_t j=0; j<i; ++j)
          {
            const size_t ptype=ig+P.GroupID[j];
            if (RecalcSwitch[ptype])
            {
              std::fill(derivs.begin(),derivs.end(),0.0);
              if (!F[ptype]->evaluateDerivatives(dist[j],derivs)) continue;
              RealType rinv(cone/dist[j]);
              PosType dr(displ[j]);
              for (int p=OffSet[ptype].first, ip=0; p<OffSet[ptype].second; ++p,++ip)
              {
                RealType dudr(rinv*derivs[ip][1]);
                RealType lap(derivs[ip][2]+lapfac*dudr);
                //RealType lap(derivs[ip][2]+(OHMMS_DIM-1.0)*dudr);
                PosType gr(dudr*dr);
                dLogPsi[p]-=derivs[ip][0];
                (*gradLogPsi[p])[i] += gr;
                (*gradLogPsi[p])[j] -= gr;
                (*lapLogPsi[p])[i] -=lap;
                (*lapLogPsi[p])[j] -=lap;
              }
            }
          }
        }
      }
      else
      {
        for (int i=0; i<d_table->size(SourceIndex); ++i)
        {
          for (int nn=d_table->M[i]; nn<d_table->M[i+1]; ++nn)
          {
            int ptype=d_table->PairID[nn];
            if (RecalcSwitch[ptype])
            {
              std::fill(derivs.begin(),derivs.end(),0.0);
              if (!F[ptype]->evaluateDerivatives(d_table->r(nn),derivs)) continue;
              int j = d_table->J[nn];
              RealType rinv(d_table->rinv(nn));
              PosType dr(d_table->dr(nn));
              for (int p=OffSet[ptype].first, ip=0; p<OffSet[ptype].second; ++p,++ip)
              {
                RealType dudr(rinv*derivs[ip][1]);
                RealType lap(derivs[ip][2]+(OHMMS_DIM-1.0)*dudr);
                PosType gr(dudr*dr);
                dLogPsi[p]-=derivs[ip][0];
                (*gradLogPsi[p])[i] += gr;
                (*gradLogPsi[p])[j] -= gr;
                (*lapLogPsi[p])[i] -=lap;
                (*lapLogPsi[p])[j] -=lap;
              }
            }
          }
        }
      }
      for (int k=0; k<myVars.size(); ++k)
      {
        int kk=myVars.where(k);
        if (kk<0)
          continue;
        if (rcsingles[k])
        {
          dlogpsi[kk]=dLogPsi[k];
          RealType dH=0.0;
          for (int iptcl=0;iptcl<NumPtcls;iptcl++)
          {
            dH += -1./(2.*mass_vec[iptcl]) * (*lapLogPsi[k])[iptcl]
              - 1./(mass_vec[iptcl]) * dot(P.G[iptcl],(*gradLogPsi[k])[iptcl]);
          }
          dhpsioverpsi[kk]=dH;
        }
        //optVars.setDeriv(p,dLogPsi[ip],-0.5*Sum(*lapLogPsi[ip])-Dot(P.G,*gradLogPsi[ip]));
      }
    }
  }

  DiffOrbitalBasePtr makeClone(ParticleSet& tqp) const
  {
    DiffTwoBodyJastrowOrbital<FT>* j2copy=new DiffTwoBodyJastrowOrbital<FT>(tqp);
    // 1. clone the map of unique pair functions "J2Unique"
    typename std::map<std::string,FT*>::const_iterator 
      it(J2Unique.begin()), it_end(J2Unique.end());
    for (;it!=it_end;it++)
    {
      FT* rf = new FT(*it->second);
      j2copy->J2Unique[it->first] = rf;
    }

    // 2. map the array of pair functions "F"
    j2copy->F.resize(F.size());
    for (int i=0;i<F.size();i++)
    { // map unique functors into the non-unique functors,
      //  also check that every pair of particles have been assigned a functor
      bool done=false;
      typename std::map<std::string,FT*>::const_iterator it(j2copy->J2Unique.begin()),it_end(j2copy->J2Unique.end());
      for (;it!=it_end;it++)
      {
        done = true;
        j2copy->F[i] = it->second;
      }
      if (!done)
      {
        APP_ABORT("Error cloning DiffTwoBodyJastrowOrbital.\n")
      }
    }
    j2copy->myVars.clear();
    j2copy->myVars.insertFrom(myVars);
    j2copy->NumVars=NumVars;
    j2copy->NumPtcls=NumPtcls;
    j2copy->NumGroups=NumGroups;
    j2copy->dLogPsi.resize(NumVars);
    j2copy->gradLogPsi.resize(NumVars,0);
    j2copy->lapLogPsi.resize(NumVars,0);
    for (int i=0; i<NumVars; ++i)
    {
      j2copy->gradLogPsi[i]=new GradVectorType(NumPtcls);
      j2copy->lapLogPsi[i]=new ValueVectorType(NumPtcls);
    }
    j2copy->OffSet=OffSet;
    return j2copy;
  }

};
}
#endif

