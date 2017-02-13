//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Yubo Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign 
//
// File created by: Yubo Yang, yyang173@illinois.edu, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCDrivers/VMC/VMCUpdateSpeciesWithDrift.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{

/// Constructor.
VMCUpdateSpeciesWithDrift::VMCUpdateSpeciesWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
    QMCHamiltonian& h, RandomGenerator_t& rg):
  QMCUpdateBase(w,psi,h,rg), tspecies(w.getSpeciesSet())
{
  // get number of species and assign defaults
  int num_species = tspecies.size();
  move_interval.assign(num_species,0); // move every step
  nrest.assign(num_species,0);         // rest counters start at 0
}

bool VMCUpdateSpeciesWithDrift::put(xmlNodePtr cur)
{
  // point a parameter set to current xml node
  ParameterSet m_param;
  m_param.put(cur);

  // initialize xml attribute reader
  OhmmsAttributeSet rattrib;
  std::string type;
  std::vector<int>  tmp_move_interval;

  // loop through child nodes
  xmlNodePtr tcur = cur->children;
  while (tcur != NULL)
  {
    std::string cname((const char*)(tcur->name));
    if (cname == "update_interval")
    {
      rattrib.add(type,"type");
      rattrib.put(tcur);
      assert(type=="Array");
      putContent(tmp_move_interval,tcur);
      assert(tmp_move_interval.size()==move_interval.size());
      copy(tmp_move_interval.begin(),tmp_move_interval.end(),move_interval.begin());
    }

    tcur = tcur->next;
  }

  // report parameters
  for (int ig=0;ig<tspecies.size();ig++){
    app_log() << " moving " << tspecies.speciesName[ig] << " every " << move_interval[ig]+1 << " steps with all-particle moves." << std::endl;
  }
  return true;
}

VMCUpdateSpeciesWithDrift::~VMCUpdateSpeciesWithDrift()
{
}

void VMCUpdateSpeciesWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
	for (; it!= it_end; ++it)
  { // get walker it and load into W
    Walker_t& thisWalker(**it);
    W.loadWalker(thisWalker,true);

    // prepare some random numbers in deltaR
    makeGaussRandomWithEngine(deltaR,RandomGen);

    for(int ig=0; ig<W.groups(); ++ig)
    { // loop over species
      if (nrest[ig]>=move_interval[ig])
      { // only move species that had enough rest

        // move species ig
        bool moved=false;
				moved = move_species_all(ig,thisWalker);

        if (moved)
        { // refresh observables
          H.auxHevaluate(W,thisWalker);
          H.saveProperty(thisWalker.getPropertyBase());
        }
        nrest[ig] = 0;
      } else { // otherwise, let them rest
        nrest[ig] += 1;
      } // if nrest > move_interval

    } // end for ig

  } // end for it (walker)
}

bool VMCUpdateSpeciesWithDrift::move_species_all(int ig,Walker_t& thisWalker)
{
  // Important data loaded in advanceWalkers:
  //   deltaR: a list of Gaussian move vectors
  //   drift:  a list of drift vectors
  //   W: local copy of thisWalker, used to store proposed configuration
  //   thisWalker: current walker about to be updated, stores reference configuration

  // prepare drift vectors ("quantum force") in drift 
  //  !! W.G must be kept up-to-date with Psi.evaluateLog(W) or else
  assignDrift(Tau,MassInvP,W.G,drift); 
  
  bool moved = false;
  int nptcl = W.last(ig) - W.first(ig);
  RealType sqrt_tau_over_mass = std::sqrt(Tau*MassInvS[ig]);
  PosType dr;        // forward move vector (i.e. R'-R)
  PosType gauss_inv; // backward gaussian vector (i.e. deltaR for reverse move)

  // calculate forward move probability
  RealType logGf=0.0; // ln T(R->R')
  for (int iat=W.first(ig); iat<W.last(ig); ++iat)
	{ // move all particles of species ig
    logGf += -0.5*dot(deltaR[iat],deltaR[iat]);
    dr = sqrt_tau_over_mass*deltaR[iat] + drift[iat];
    W.makeMoveAndCheck(iat,dr); // move and check PBC
    W.acceptMove(iat); // update distance tables and SK
	} //iat

  // !! must be called before assignDrift for backward move
  RealType logpsi  = Psi.evaluateLog(W); // !! side effects: updates W.G and W.L, also potentialy W.SK
  // calculate reverse move probability
  //  !! W.G must be kept up-to-date with Psi.evaluateLog(W) or else
  assignDrift(Tau,MassInvP,W.G,drift); // assign backward drifts
  RealType logGb=0.0; // ln T(R'->R)
  for (int iat=W.first(ig); iat<W.last(ig); ++iat)
  {
    gauss_inv = (thisWalker.R[iat] - W.R[iat] - drift[iat])/sqrt_tau_over_mass;
    logGb += -0.5*dot(gauss_inv,gauss_inv);
  }

  // accept/reject
  RealType logpsi0 = thisWalker.Properties(LOGPSI);
  RealType prob = std::exp(logGb-logGf+2.0*(logpsi-logpsi0));
  if (RandomGen() > prob)
  {
    nReject += nptcl; // expand to # of single particle moves
  } else {
    moved = true;
    nAccept += nptcl; // expand to # of single particle moves
  }

  if (moved)
  { // transfer walker changes
    W.saveWalker(thisWalker);
    RealType phase = Psi.getPhase();
    EstimatorRealType eloc = H.evaluate(W); // !! side effects: updates W.PropertyList[LOCALENERGY] and P.PropertyList[LOCALPOTENTIAL]
    thisWalker.resetProperty(logpsi,phase,eloc); // !! side effect: thisWalker.Age=0
  } else {
    H.rejectedMove(W,thisWalker);
  }

  return moved;
}
}

/***************************************************************************
 * $RCSfile: VMCUpdateSpeciesWithDrift.cpp $   $Author: yyang173 $
 * $Revision: 7223 $   $Date: 2016/11/14 14:29:40 $
 * $Id: VMCUpdateSpecies.h,v 1.5 2016/11/14 14:29:40 yyang173 Exp $
 ***************************************************************************/
