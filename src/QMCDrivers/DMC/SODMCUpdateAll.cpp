#include "SODMCUpdateAll.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "QMCDrivers/DriftOperators.h"
#include "QMCDrivers/WalkerProperties.h"

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

/// Constructor.
SODMCUpdateAll::SODMCUpdateAll(MCWalkerConfiguration& w,
                                                     TrialWaveFunction& psi,
                                                     QMCHamiltonian& h,
                                                     RandomGenerator& rg)
    : QMCUpdateBase(w, psi, h, rg)
{
  UpdatePbyP = false;
}

/// destructor
SODMCUpdateAll::~SODMCUpdateAll() {}

void SODMCUpdateAll::advanceWalker(Walker_t& thisWalker, bool recompute)
{
  W.loadWalker(thisWalker, false);
  //create a 3N-Dimensional Gaussian with variance=1
  RealType nodecorr = setScaledDriftPbyPandNodeCorr(Tau, MassInvP, W.G, drift);
  //RealType nodecorr = setScaledDriftPbyPandNodeCorr(m_tauovermass,W.G,drift);
  makeGaussRandomWithEngine(deltaR, RandomGen);
  if (ndim < 3)
  {
    for (int iat = 0; iat < deltaR.size(); ++iat)
      deltaR[iat][2] = 0;
  }
  makeGaussRandomWithEngine(deltaS, RandomGen);
  //save old local energy
  RealType eold        = thisWalker.Properties(WP::LOCALENERGY);
  RealType enew        = eold;
  bool accepted        = false;
  RealType rr_accepted = 0.0;
  RealType rr_proposed = 0.0;
  RealType logpsi;

  RealType invSqrtSpinMass = 1.0 / std::sqrt(spinMass);

  if (W.makeMoveAllParticlesWithDrift(thisWalker, drift, deltaR, SqrtTauOverMass))
  {
    for (int iat = 0; iat < deltaS.size(); ++iat)
      W.spins[iat] = thisWalker.spins[iat] + invSqrtSpinMass * SqrtTauOverMass[iat] * deltaS[iat];
    //evaluate the new wave function
    logpsi = Psi.evaluateLog(W);
    //fixed node
    if (!branchEngine->phaseChanged(Psi.getPhaseDiff()))
    {
      RealType logGf = -0.5 * Dot(deltaR, deltaR);
      nodecorr       = setScaledDriftPbyPandNodeCorr(Tau, MassInvP, W.G, drift);
      deltaR         = thisWalker.R - W.R - drift;
      RealType logGb = logBackwardGF(deltaR);
      //RealType logGb = -m_oneover2tau*Dot(deltaR,deltaR);
      RealType prob = std::min(std::exp(logGb - logGf + 2.0 * (logpsi - thisWalker.Properties(WP::LOGPSI))), 1.0);
      //calculate rr_proposed here
      deltaR      = W.R - thisWalker.R;
      rr_proposed = Dot(deltaR, deltaR);
      if (RandomGen() <= prob)
      {
        accepted    = true;
        rr_accepted = rr_proposed;
      }
    }
  }

  // recompute Psi if the move is rejected
  if (!accepted)
  {
    W.loadWalker(thisWalker, false);
    W.update();
    logpsi = Psi.evaluateLog(W);
  }

  // evaluate Hamiltonian
  enew = H.evaluateWithToperator(W);
  H.auxHevaluate(W, thisWalker);
  H.saveProperty(thisWalker.getPropertyBase());

  // operate on thisWalker.
  if (accepted)
  {
    W.saveWalker(thisWalker);
    thisWalker.resetProperty(logpsi, Psi.getPhase(), enew, rr_accepted, rr_proposed, nodecorr);
  }
  else
  {
    thisWalker.Age++;
    thisWalker.Properties(WP::R2ACCEPTED) = 0.0;
    thisWalker.Properties(WP::R2PROPOSED) = rr_proposed;
  }

  const int NonLocalMoveAcceptedTemp = H.makeNonLocalMoves(W);
  if (NonLocalMoveAcceptedTemp > 0)
  {
    W.saveWalker(thisWalker);
    thisWalker.resetProperty(Psi.getLogPsi(), Psi.getPhase(), enew);
    // debugging lines
    //logpsi = Psi.getLogPsi();
    //W.update(true);
    //RealType logpsi2 = Psi.evaluateLog(W);
    //if(logpsi!=logpsi2) std::cout << " logpsi " << logpsi << " logpsi2 " << logpsi2
    //                              << " diff " << logpsi2-logpsi << std::endl;

    NonLocalMoveAccepted += NonLocalMoveAcceptedTemp;
  }

  thisWalker.Weight *= branchEngine->branchWeight(enew, eold);
  //branchEngine->accumulate(eold,1);
  if (accepted)
    ++nAccept;
  else
    ++nReject;

  setMultiplicity(thisWalker);
}
} // namespace qmcplusplus
