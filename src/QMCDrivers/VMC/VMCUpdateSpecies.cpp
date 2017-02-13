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

#include "QMCDrivers/VMC/VMCUpdateSpecies.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{

/// Constructor.
VMCUpdateSpecies::VMCUpdateSpecies(MCWalkerConfiguration& w, TrialWaveFunction& psi,
    QMCHamiltonian& h, RandomGenerator_t& rg):
  QMCUpdateBase(w,psi,h,rg), tspecies(w.getSpeciesSet())
{
  //UpdatePbyP=true; // this allows SqrtTauOverMass to be filled in QMCUpdateBase::resetRun
  // took out flag in QMCUpdateBase

  // get number of species and assign defaults
  int num_species = tspecies.size();
  move_interval.assign(num_species,0); // move every step
  move_pbyp.assign(num_species,false); // move with pbyp=false, eaiser to get right (no w_buffer)
  nrest.assign(num_species,0);         // rest counters start at 0
}

VMCUpdateSpecies::~VMCUpdateSpecies()
{
}

bool VMCUpdateSpecies::put(xmlNodePtr cur)
{
  // point a parameter set to current xml node
  ParameterSet m_param;
  m_param.put(cur);

  // initialize xml attribute reader
  OhmmsAttributeSet rattrib;
  std::string type;
  std::vector<int>  tmp_move_interval;
  std::vector<bool> tmp_move_pbyp;

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
    } else if (cname == "species_pbyp") 
    {
      rattrib.add(type,"type");
      rattrib.put(tcur);
      assert(type=="Array");
      putContent(tmp_move_pbyp,tcur);
      assert(tmp_move_pbyp.size()==move_pbyp.size());
      copy(tmp_move_pbyp.begin(),tmp_move_pbyp.end(),move_pbyp.begin());
    }

    tcur = tcur->next;
  }

  // report parameters
  for (int ig=0;ig<tspecies.size();ig++){
    if (move_pbyp[ig])
    {
      APP_ABORT("pbyp currently not properly implemented");
    }
    app_log() << " moving " << tspecies.speciesName[ig] << " every " << move_interval[ig]+1 << " steps with pbyp=" << move_pbyp[ig] << std::endl;
  }
  return true;
}

void VMCUpdateSpecies::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  for (; it!= it_end; ++it)
  { // get walker it and load into W
    Walker_t& thisWalker(**it);
    W.loadWalker(thisWalker,true);

    // YY: not sure why w_buffer is needed, but it seems to be important for pbyp
    //  will try to keep w_buffer up-to-date with Psi.updateBuffer whenever 
    //  moves are accepted.
    // DataSet is initialized by QMCUpdateBase::initWalkersForPbyP
    //  , initWalkersForPbyP also asked Psi to track thisWalker.DataSet
    Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
    Psi.copyFromBuffer(W,w_buffer);

    // prepare some random numbers in deltaR
    makeGaussRandomWithEngine(deltaR,RandomGen);

    for(int ig=0; ig<W.groups(); ++ig) 
    { // loop over species
      if (nrest[ig]>=move_interval[ig])
      { // only move species that had enough rest

        // move species ig with selected method
        //  thisWalker and w_buffer will be updated
        bool moved=false;
        if (move_pbyp[ig])
        {
          moved = move_species_pbyp(ig,thisWalker,w_buffer);
        } else {
          moved = move_species_all(ig,thisWalker,w_buffer);
        }

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

bool VMCUpdateSpecies::move_species_pbyp(int ig,Walker_t& thisWalker, Walker_t::Buffer_t& w_buffer)
{
  // W: MCWalkerConfiguration loaded in advanceWalkers
  // deltaR: a list of random numbers loaded in advanceWalkers
  // thisWalker: current walker about to be updated, contains reference data

  bool moved_one = false;
  RealType sqrttau = std::sqrt(Tau*MassInvS[ig]);
  for (int iat=W.first(ig); iat<W.last(ig); ++iat)
	{ // make single particle move for each ptcl in species ig
		mPosType dr = sqrttau*deltaR[iat];
		if (W.makeMoveAndCheck(iat,dr))
		{ // if move was legal
			RealType ratio = Psi.ratio(W,iat);
			RealType prob  = ratio*ratio;
			if (RandomGen() < prob)
			{ // accept move
        moved_one = true;
				++nAccept;
				W.acceptMove(iat);
				Psi.acceptMove(W,iat);
			} else { // reject move
				++nReject;
				W.rejectMove(iat);
				Psi.rejectMove(iat);
			}
		} else { // reject illegal moves
			++nReject;
    }
	} //iat

  // transfer walker changes
  W.saveWalker(thisWalker);
  // update w_buffer (false: not from scratch)
  RealType logpsi = Psi.updateBuffer(W,w_buffer,false);

  if (moved_one)
  {
		EstimatorRealType eloc=H.evaluate(W);
		thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
  }

  return moved_one;
}

bool VMCUpdateSpecies::move_species_all(int ig,Walker_t& thisWalker, Walker_t::Buffer_t& w_buffer)
{
  // W: MCWalkerConfiguration loaded in advanceWalkers
  // deltaR: a list of random numbers loaded in advanceWalkers
  // thisWalker: current walker about to be updated, contains reference data

  int nptcl = W.last(ig) - W.first(ig);
  RealType sqrttau = std::sqrt(Tau*MassInvS[ig]);
  for (int iat=W.first(ig); iat<W.last(ig); ++iat)
	{ // move all particles of species ig
		mPosType dr = sqrttau*deltaR[iat];
    W.makeMoveAndCheck(iat,dr); // move and check PBC
    W.acceptMove(iat); // update distance tables and SK
    //Psi.acceptMove(W,iat); // do NOT update the wave function value, not actually accepting these moves yet
	} //iat

  bool moved = false;
  RealType logpsi(Psi.evaluateLog(W));
  RealType prob = std::exp(2.0*(logpsi-thisWalker.Properties(LOGPSI)));
  if (RandomGen() > prob)
  {
    //thisWalker.Age++; // YY: necessary?
    nReject += nptcl; // expand to # of single particle moves
    H.rejectedMove(W,thisWalker);
  } else {
    moved = true;
    nAccept += nptcl; // expand to # of single particle moves
  }

  if (moved)
	{ 
    // transfer walker changes
    W.saveWalker(thisWalker);

    // MUST NOT update buffer if not using pbyp ?!
    // update w_buffer from scratch
    //RealType logpsi1 = Psi.updateBuffer(W,w_buffer,true);
    // logpsi1 == logpsi ?
    
    // MUST NOT reset walker every step if not using pbyp ?!
    EstimatorRealType eloc=H.evaluate(W);
    thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
	}

  return moved;
}

}

/***************************************************************************
 * $RCSfile: VMCUpdateSpecies.h,v $   $Author: yyang173 $
 * $Revision: 7223 $   $Date: 2016/11/14 14:29:40 $
 * $Id: VMCUpdateSpecies.h,v 1.5 2016/11/14 14:29:40 yyang173 Exp $
 ***************************************************************************/
