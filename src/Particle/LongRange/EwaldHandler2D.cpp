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

#include "Particle/LongRange/EwaldHandler2D.h"

namespace qmcplusplus
{

EwaldHandler2D::EwaldHandler2D(ParticleSet& ref, mRealType kc_in)
  : LRHandlerBase(kc_in)
{
  if (ref.Lattice.ndim != 2)
    throw std::runtime_error("2D Ewald requires 2D Lattice");
  LR_rc = ref.Lattice.LR_rc; // CoulombPBC needs get_rc() to createSpline4RbyVs
  LR_kc = ref.Lattice.LR_kc; // get_kc() is used in QMCFiniteSize
  alpha = std::sqrt(LR_kc/2.0/LR_rc);
  area = ref.Lattice.Volume/ref.Lattice.R(2,2);
  // report
  app_log() << "    alpha = " << alpha << " area = " << area << std::endl;
  fillFk(ref.SK->KLists);
  fillZheights(ref);
}

void EwaldHandler2D::fillFk(const KContainer& KList)
{
  const mRealType knorm = M_PI / area;
  const mRealType kalpha = 1.0 / (2.0*alpha);
  mRealType kmag, uk;

  Fk.resize(KList.kpts_cart.size());
  MaxKshell = KList.kshell.size() - 1;
  Fk_symm.resize(MaxKshell);
  kmags.resize(MaxKshell);

  for (int ks = 0, ki = 0; ks < Fk_symm.size(); ks++)
  {
    kmag = std::sqrt(KList.ksq[ki]);
    kmags[ks] = kmag;
    uk = knorm/kmag;
    Fk_symm[ks] = uk;
    while (ki < KList.kshell[ks + 1] && ki < Fk.size())
      Fk[ki++] = uk;
  }
}

void EwaldHandler2D::fillZheights(const ParticleSet& P)
{
  const SpeciesSet& tspecies(P.getSpeciesSet());
  const int nspec = tspecies.TotalNum;
  if (nspec > 4)
    throw std::runtime_error("too many species");
  mRealType zij0, zij;
  for (int ispec=0; ispec<nspec; ispec++)
  {
    for (int jspec=0; jspec<nspec; jspec++)
    {
      zij0 = std::abs(P.R[P.first(ispec)][2]-P.R[P.first(jspec)][2]);
      for (int i=P.first(ispec); i<P.last(ispec); i++)
      {
        for (int j=P.first(jspec); j<P.last(jspec); j++)
        {
          zij = std::abs(P.R[i][2]-P.R[j][2]);
          if (std::abs(zij-zij0)>1e-12)
            throw std::runtime_error("species not in xy plane");
        }
      }
      zheights[ispec, jspec] = zij0;
    }
  }
}

EwaldHandler2D::mRealType EwaldHandler2D::slab_func(mRealType z, mRealType k) const
{
  mRealType term1, term2;
  term1 = std::exp(slab_logf( z, k));
  term2 = std::exp(slab_logf(-z, k));
  return term1+term2;
}

EwaldHandler2D::mRealType EwaldHandler2D::slab_logf(mRealType z, mRealType k) const
{
  mRealType z1 = alpha*z;
  mRealType k1 = k/(2*alpha);
  mRealType log_term = k*z + std::log(erfc(z1+k1));
  return log_term;
}

EwaldHandler2D::mRealType EwaldHandler2D::slab_vsr_k0(mRealType z) const
{
  mRealType z1 = alpha*z;
  mRealType term1 = std::exp(-z1*z1) * 2*std::sqrt(M_PI)/alpha;
  mRealType term2 = 2*M_PI*z*erfc(z1);
  return term1-term2;
}


EwaldHandler2D::mRealType EwaldHandler2D::evaluateLayers(
  const std::vector<int>& kshell,
  const pRealType* restrict rk1_r,
  const pRealType* restrict rk1_i,
  const pRealType* restrict rk2_r,
  const pRealType* restrict rk2_i,
  const int ispec, const int jspec
) const
{
  mRealType z = zheights[ispec, jspec];
  mRealType vk = 0.0;
  mRealType uk = 0.0;
  for (int ks = 0, ki = 0; ks < MaxKshell; ks++)
  {
    mRealType u = 0;
    for (; ki < kshell[ks + 1]; ki++)
    {
      u = ((*rk1_r++) * (*rk2_r++) + (*rk1_i++) * (*rk2_i++));
      uk = Fk[ki]*slab_func(z, kmags[ks]);
      vk += u * uk;
    }
  }
  return vk;
}

EwaldHandler2D::mRealType EwaldHandler2D::sumMadelung(
  const std::vector<int>& kshell
) const
{
  mRealType vk = 0.0, uk;
  for (int ks = 0; ks < MaxKshell; ks++)
  {
    const int nterm = kshell[ks + 1] - kshell[ks];
    uk = Fk_symm[ks]*slab_func(0, kmags[ks]);
    vk += nterm*uk;
  }
  return vk;
}

EwaldHandler2D::mRealType EwaldHandler2D::evaluateBackground(
  const ParticleSet& P,
  const int ispec,
  const int jspec) const
{
  mRealType z = zheights[ispec, jspec];
  return slab_vsr_k0(z)/area;
}

} // qmcplusplus
