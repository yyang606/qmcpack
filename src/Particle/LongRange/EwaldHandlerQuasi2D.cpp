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

#include "Particle/LongRange/EwaldHandlerQuasi2D.h"

namespace qmcplusplus
{

EwaldHandlerQuasi2D::EwaldHandlerQuasi2D(ParticleSet& ref, mRealType kc_in)
  : LRHandlerBase(kc_in)
{
  if (ref.Lattice.ndim != 2)
    throw std::runtime_error("Quasi2D Ewald requires 2D Lattice");
  LR_rc = ref.Lattice.LR_rc; // CoulombPBC needs get_rc() to createSpline4RbyVs
  LR_kc = ref.Lattice.LR_kc; // get_kc() is used in QMCFiniteSize
  alpha = std::sqrt(LR_kc/2.0/LR_rc);
  area = ref.Lattice.Volume/ref.Lattice.R(2,2);
  // report
  app_log() << "    alpha = " << alpha << " area = " << area << std::endl;
  fillFk(ref.SK->KLists);
}

void EwaldHandlerQuasi2D::fillFk(const KContainer& KList)
{
  const mRealType knorm = M_PI / area;
  mRealType kmag, uk;

  Fk.resize(KList.kpts_cart.size());
  MaxKshell = KList.kshell.size() - 1;
  Fk_symm.resize(MaxKshell);

  kMag.resize(MaxKshell);
  for (int ks = 0, ki = 0; ks < Fk_symm.size(); ks++)
  {
    kmag = std::sqrt(KList.ksq[ki]);
    kMag[ks] = kmag; // store k magnitutes
    uk = knorm/kmag;
    Fk_symm[ks] = uk;
    while (ki < KList.kshell[ks + 1] && ki < Fk.size())
      Fk[ki++] = uk;
  }
}

EwaldHandlerQuasi2D::mRealType EwaldHandlerQuasi2D::evaluate_slab(
  pRealType z, const std::vector<int>& kshell,
  const pComplexType* restrict eikr_i,
  const pComplexType* restrict eikr_j)
{
  mRealType zmag = std::abs(z);
  mRealType vk = -slab_vsr_k0(zmag)/area;
  for (int ks = 0, ki = 0; ks < MaxKshell; ks++)
  {
    mRealType u = 0;
    for (; ki < kshell[ks + 1]; ki++, eikr_i++, eikr_j++)
      u += ((*eikr_i).real() * (*eikr_j).real() + (*eikr_i).imag() * (*eikr_j).imag());
    vk += u * Fk_symm[ks] * slab_func(zmag, kMag[ks]);
  }
  return vk;
}

EwaldHandlerQuasi2D::mRealType EwaldHandlerQuasi2D::evaluate_slab(
  pRealType z, const std::vector<int>& kshell,
  const pRealType* restrict rk1_r,
  const pRealType* restrict rk1_i,
  const pRealType* restrict rk2_r,
  const pRealType* restrict rk2_i) const
{
  mRealType zmag = std::abs(z);
  mRealType vk = -slab_vsr_k0(zmag)/area;
  for (int ks = 0, ki = 0; ks < MaxKshell; ks++)
  {
    mRealType u = 0;
    for (; ki < kshell[ks + 1]; ki++)
      u += (*rk1_r++) * (*rk2_r++) + (*rk1_i++) * (*rk2_i++);
    vk += u * Fk_symm[ks] * slab_func(zmag, kMag[ks]);
  }
  return vk;
}

EwaldHandlerQuasi2D::mRealType EwaldHandlerQuasi2D::slab_func(mRealType z, mRealType k) const
{
  mRealType term1, term2;
  term1 = std::exp(slab_logf( z, k));
  term2 = std::exp(slab_logf(-z, k));
  return term1+term2;
}

EwaldHandlerQuasi2D::mRealType EwaldHandlerQuasi2D::slab_logf(mRealType z, mRealType k) const
{
  mRealType z1 = alpha*z;
  mRealType k1 = k/(2*alpha);
  mRealType log_term = k*z + std::log(erfc(z1+k1));
  return log_term;
}

EwaldHandlerQuasi2D::mRealType EwaldHandlerQuasi2D::slab_vsr_k0(mRealType z) const
{
  mRealType z1 = alpha*z;
  mRealType term1 = std::exp(-z1*z1) * 2*std::sqrt(M_PI)/alpha;
  mRealType term2 = 2*M_PI*z*erfc(z1);
  return term1-term2;
}

} // qmcplusplus
