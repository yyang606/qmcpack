//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#ifndef QMCPLUSPLUS_OLD2D_HANDLER_H
#define QMCPLUSPLUS_OLD2D_HANDLER_H

#include "LongRange/LRHandlerBase.h"

namespace qmcplusplus
{
/* LR breakup for the standard Ewald method
 *
 * Quasi-2D Ewald method : J. Phys.: Condens. Matter 16, 891 (2004)
 * http://iopscience.iop.org/0953-8984/16/6/017/
 * Note that \f$ \simga \rightarrow 1/\sigma\f$
 * It is possible to use 3D Ewald but for the bulk system, the optimal breakup method
 * is used.
 */
class TwoDEwaldHandler : public LRHandlerBase
{
public:
  /// Related to the Gaussian width: \f$ v_l = v(r)erf(\sigma r)\f$
  mRealType Sigma;
  mRealType Volume;
  ///store |k|
  std::vector<mRealType> kMag;
  /// Constructor
  TwoDEwaldHandler(ParticleSet& ref, mRealType kc_in = -1.0) : LRHandlerBase(kc_in)
  {
    Sigma = LR_kc = ref.Lattice.LR_kc;
  }

  /** "copy" constructor
   * @param aLR LRHandlerTemp
   * @param ref Particleset
   *
   * Copy the content of aLR
   * References to ParticleSet or ParticleLayoutout_t are not copied.
   */
  TwoDEwaldHandler(const TwoDEwaldHandler& aLR, ParticleSet& ref);

  LRHandlerBase* makeClone(ParticleSet& ref) override { return new TwoDEwaldHandler(*this, ref); }

  void initBreakup(ParticleSet& ref) override;

  void Breakup(ParticleSet& ref, mRealType rs_in) override { initBreakup(ref); }

  void resetTargetParticleSet(ParticleSet& ref) override {}

  inline mRealType evaluate(mRealType r, mRealType rinv) override { return erfc(r * Sigma) * rinv; }

  /** evaluate the contribution from the long-range part for for spline
   */
  inline mRealType evaluateLR(mRealType r) override { return -erf(r * Sigma) / r; }

  inline mRealType evaluate_vlr_k(mRealType k) override { return 2.0*std::sqrt(M_PI)/k * erfc(k/(2*Sigma)) / Volume; }

  inline mRealType evaluateSR_k0() override { return 2.0 * std::sqrt(M_PI) / (Sigma * Volume); }

  inline mRealType evaluateLR_r0() override { return 2.0 * Sigma / std::sqrt(M_PI); }

  /**  evaluate the first derivative of the short range part at r
   *
   * @param r  radius
   * @param rinv 1/r
   */
  inline mRealType srDf(mRealType r, mRealType rinv) override { return 0.0; }

  void fillFk(KContainer& KList);
};
} // namespace qmcplusplus
#endif
