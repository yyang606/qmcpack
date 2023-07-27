//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_EWALDSRSCREEN2D_H
#define QMCPLUSPLUS_EWALDSRSCREEN2D_H

#include "Particle/LongRange/LRHandlerBase.h"

namespace qmcplusplus
{
class Screen2D : public LRHandlerBase
{
public:
  Screen2D(ParticleSet& ref, mRealType kc_in=-1.0);

  // copy constructor
  LRHandlerBase* makeClone(ParticleSet& ref) const override { return new Screen2D(*this); }

  // short-range part
  mRealType evaluate(mRealType r, mRealType rinv) const override;
  inline mRealType evaluateLR_r0() const override {return 0.0;}

  // long-range part
  inline mRealType evaluateLR(mRealType r) const override {return 0.0;}
  mRealType evaluateSR_k0() const override;
  void fillFk(const KContainer& KList);

  // begin required overrides
  inline mRealType srDf(mRealType r, mRealType rinv) const override
  {
    throw std::runtime_error("2D Ewald sr derivative not implemented");
  }
  inline mRealType lrDf(mRealType r) const override
  {
    throw std::runtime_error("2D Ewald lr derivative not implemented");
  }
  virtual mRealType evaluate_vlr_k(mRealType k) const override
  {
    throw std::runtime_error("2D Ewald vlr_k not implemented");
  }
  void initBreakup(ParticleSet& ref) override {}
  void Breakup(ParticleSet& ref, mRealType rs_in) override { initBreakup(ref); }
  void resetTargetParticleSet(ParticleSet& ref) override {}
  // overrides end
private:
  mRealType area;
  mRealType dgate;
  int mimg;
};
} // qmcplusplus
#endif
