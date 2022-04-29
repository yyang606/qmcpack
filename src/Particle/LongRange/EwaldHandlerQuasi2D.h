//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_EWALD_QUASI2D_HANDLER_H
#define QMCPLUSPLUS_EWALD_QUASI2D_HANDLER_H

#include "Particle/LongRange/LRHandlerBase.h"

namespace qmcplusplus
{
/* LR breakup for the standard Ewald method in quasi-2D (slab geometry)
 */
class EwaldHandlerQuasi2D : public LRHandlerBase
{
public:
  EwaldHandlerQuasi2D(ParticleSet& ref, mRealType kc_in=-1.0);

  // copy constructor
  LRHandlerBase* makeClone(ParticleSet& ref) override { return new EwaldHandlerQuasi2D(*this); }

  // short-range part
  inline mRealType evaluate(mRealType r, mRealType rinv) override { return erfc(alpha*r) * rinv; }
  inline mRealType evaluateLR_r0() override { return 2.0 * alpha / std::sqrt(M_PI); }

  // long-range part
  inline mRealType evaluateLR(mRealType r) override { return erf(alpha*r) / r; }
  inline mRealType evaluateSR_k0() override { return slab_vsr_k0(0)/area; }
  void fillFk(const KContainer& KList);

  // z-dependent long-range part
  mRealType evaluate_slab(pRealType z,
                          const std::vector<int>& kshell,
                          const pComplexType* restrict eikr_i,
                          const pComplexType* restrict eikr_j) override;

  // begin required overrides
  inline mRealType srDf(mRealType r, mRealType rinv) override
  {
    throw std::runtime_error("Quasi2D Ewald sr derivative not implemented");
  }
  inline mRealType lrDf(mRealType r) override
  {
    throw std::runtime_error("Quasi2D Ewald lr derivative not implemented");
  }
  virtual mRealType evaluate_vlr_k(mRealType k) override
  {
    throw std::runtime_error("Quasi2D Ewald vlr_k not implemented");
  }
  void initBreakup(ParticleSet& ref) override {}
  void Breakup(ParticleSet& ref, mRealType rs_in) override { initBreakup(ref); }
  void resetTargetParticleSet(ParticleSet& ref) override {}
  // overrides end
private:
  mRealType alpha;
  mRealType area;
  ///store |k|
  std::vector<mRealType> kMag;
  mRealType slab_func(mRealType z, mRealType k) const;
  mRealType slab_logf(mRealType z, mRealType k) const;
  mRealType slab_vsr_k0(mRealType z) const;
};
} // qmcplusplus
#endif
