//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCDrivers/DMC/DMCUpdatePbyPVMC.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif
//#define TEST_INNERBRANCH
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"


namespace qmcplusplus
{
//TimerNameList_t<DMCTimers> DMCTimerNames = {{DMC_buffer, "DMCUpdatePbyP::Buffer"},
//                                            {DMC_movePbyP, "DMCUpdatePbyP::movePbyP"},
//                                            {DMC_hamiltonian, "DMCUpdatePbyP::Hamiltonian"},
//                                            {DMC_collectables, "DMCUpdatePbyP::Collectables"},
//                                            {DMC_tmoves, "DMCUpdatePbyP::Tmoves"}};


/// Constructor.
DMCUpdatePbyPVMC::DMCUpdatePbyPVMC(MCWalkerConfiguration& w,
                                                               TrialWaveFunction& psi,
                                                               QMCHamiltonian& h,
                                                               RandomGenerator_t& rg)
    : QMCUpdateBase(w, psi, h, rg)
{
  setup_timers(myTimers, DMCTimerNames, timer_level_medium);
}

/// destructor
DMCUpdatePbyPVMC::~DMCUpdatePbyPVMC() {}

void DMCUpdatePbyPVMC::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  bool accept_reject=false;

  myTimers[DMC_buffer]->start();
  Walker_t::WFBuffer_t& w_buffer(thisWalker.DataSet);
  W.loadWalker(thisWalker, true);
  Psi.copyFromBuffer(W, w_buffer);
  myTimers[DMC_buffer]->stop();
  //create a 3N-Dimensional Gaussian with variance=1
  makeGaussRandomWithEngine(deltaR, RandomGen);
  int nAcceptTemp(0);
  int nRejectTemp(0);
  //copy the old energy and scale factor of drift
  EstimatorRealType eold(thisWalker.Properties(LOCALENERGY));
  EstimatorRealType enew(eold);
  RealType rr_proposed = 0.0;
  RealType rr_accepted = 0.0;
  RealType gf_acc      = 1.0;
  myTimers[DMC_movePbyP]->start();
  for (int ig = 0; ig < W.groups(); ++ig) //loop over species
  {
    RealType tauovermass = Tau * MassInvS[ig];
    RealType oneover2tau = 0.5 / (tauovermass);
    RealType sqrttau     = std::sqrt(tauovermass);
    for (int iat = W.first(ig); iat < W.last(ig); ++iat)
    {
      W.setActive(iat);
      //get the displacement
      GradType grad_iat = Psi.evalGrad(W, iat);
      mPosType dr;
      // VMC directly from sampling, use bare unscaled drift
      if(accept_reject)
        getScaledDrift(tauovermass, grad_iat, dr);
      else
        getUnscaledDrift(tauovermass, grad_iat, dr);
      dr += sqrttau * deltaR[iat];
      //RealType rr=dot(dr,dr);
      RealType rr = tauovermass * dot(deltaR[iat], deltaR[iat]);
      rr_proposed += rr;
      // VMC directly from sampling, no accept/reject
      if(accept_reject)
        if (rr > m_r2max)
        {
          ++nRejectTemp;
          continue;
        }
      
      // VMC directly from sampling, force makeMoveAndCheck to always take action and return true
      if(!accept_reject)
        W.UseBoundBox=false;
      if (!W.makeMoveAndCheck(iat, dr))
        continue;
      RealType ratio = Psi.ratioGrad(W, iat, grad_iat);

      // VMC directly from sampling, no accept/reject (i.e. always accept)
      if(!accept_reject)
      {
        ++nAcceptTemp;
        Psi.acceptMove(W, iat);
        W.acceptMove(iat);
        rr_accepted += rr;
        gf_acc *= 1.0; //acceptance ratio is 1.0
      }
      else
      {
        //node is crossed reject the move
        //if(Psi.getPhase() > std::numeric_limits<RealType>::epsilon())
        //if(branchEngine->phaseChanged(Psi.getPhase(),thisWalker.Properties(SIGN)))
        if (branchEngine->phaseChanged(Psi.getPhaseDiff()))
        {
          ++nRejectTemp;
          ++nNodeCrossing;
          W.rejectMove(iat);
          Psi.rejectMove(iat);
        }
        else
        {
          EstimatorRealType logGf = -0.5 * dot(deltaR[iat], deltaR[iat]);
          //Use the force of the particle iat
          //RealType scale=getDriftScale(m_tauovermass,grad_iat);
          //dr = W.R[iat]-W.activePos-scale*real(grad_iat);
          getScaledDrift(tauovermass, grad_iat, dr);
          dr                      = W.R[iat] - W.activePos - dr;
          EstimatorRealType logGb = -oneover2tau * dot(dr, dr);
          RealType prob           = ratio * ratio * std::exp(logGb - logGf);
          if (RandomGen() < prob)
          {
            ++nAcceptTemp;
            Psi.acceptMove(W, iat);
            W.acceptMove(iat);
            rr_accepted += rr;
            gf_acc *= prob; //accumulate the ratio
          }
          else
          {
            ++nRejectTemp;
            W.rejectMove(iat);
            Psi.rejectMove(iat);
          }
        }
      }
    }
  }
  Psi.completeUpdates();
  W.donePbyP();
  myTimers[DMC_movePbyP]->stop();

  if (nAcceptTemp > 0)
  {
    //need to overwrite the walker properties
    myTimers[DMC_buffer]->start();
    thisWalker.Age  = 0;
    RealType logpsi = Psi.updateBuffer(W, w_buffer, recompute);
    W.saveWalker(thisWalker);
    myTimers[DMC_buffer]->stop();
    myTimers[DMC_hamiltonian]->start();
    enew = H.evaluateWithToperator(W);
    myTimers[DMC_hamiltonian]->stop();
    thisWalker.resetProperty(logpsi, Psi.getPhase(), enew, rr_accepted, rr_proposed, 1.0);
    //thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
    thisWalker.Weight *= 1.0;
    myTimers[DMC_collectables]->start();
    H.auxHevaluate(W, thisWalker);
    H.saveProperty(thisWalker.getPropertyBase());
    myTimers[DMC_collectables]->stop();
  }
  else
  {
    //all moves are rejected: does not happen normally with reasonable wavefunctions
    thisWalker.Age++;
    thisWalker.Properties(R2ACCEPTED) = 0.0;
    //weight is set to 0 for traces
    // consistent w/ no evaluate/auxHevaluate
    RealType wtmp     = thisWalker.Weight;
    thisWalker.Weight = 0.0;
    H.rejectedMove(W, thisWalker);
    thisWalker.Weight = wtmp;
    ++nAllRejected;
    enew   = eold; //copy back old energy
    gf_acc = 1.0;
    //thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
    thisWalker.Weight *= 1.0;
  }
#if !defined(REMOVE_TRACEMANAGER)
  Traces->buffer_sample(W.current_step);
#endif
  myTimers[DMC_tmoves]->start();
  const int NonLocalMoveAcceptedTemp = H.makeNonLocalMoves(W);
  if (NonLocalMoveAcceptedTemp > 0)
  {
    RealType logpsi = Psi.updateBuffer(W, w_buffer, false);
    // debugging lines
    //W.update(true);
    //RealType logpsi2 = Psi.evaluateLog(W);
    //if(logpsi!=logpsi2) std::cout << " logpsi " << logpsi << " logps2i " << logpsi2 << " diff " << logpsi2-logpsi << std::endl;
    W.saveWalker(thisWalker);
    NonLocalMoveAccepted += NonLocalMoveAcceptedTemp;
  }
  myTimers[DMC_tmoves]->stop();
  nAccept += nAcceptTemp;
  nReject += nRejectTemp;

  setMultiplicity(thisWalker);
}

} // namespace qmcplusplus
