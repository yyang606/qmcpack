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

#include "Particle/LongRange/EwaldScreen2D.h"

namespace qmcplusplus
{

EwaldScreen2D::EwaldScreen2D(ParticleSet& ref, mRealType kc_in)
  : LRHandlerBase(kc_in), dgate(ref.getLattice().dgate)
{
  if (ref.getLattice().ndim != 2)
    throw std::runtime_error("2D Ewald requires 2D Lattice");
  if (dgate < 0)
    throw std::runtime_error("screened Ewald requires distance_to_gate");
  LR_rc = ref.getLattice().LR_rc; // CoulombPBC needs get_rc() to createSpline4RbyVs
  LR_kc = ref.getLattice().LR_kc; // get_kc() is used in QMCFiniteSize
  alpha = std::sqrt(LR_kc/2.0/LR_rc);
  area = ref.getLattice().Volume/ref.getLattice().R(2,2);
  // report
  app_log() << "    alpha = " << alpha << " area = " << area << std::endl;
  app_log() << "    dgate = " << dgate << std::endl;
  fillFk(ref.getSimulationCell().getKLists());
}

void EwaldScreen2D::fillFk(const KContainer& KList)
{
  const mRealType knorm = 2.0*M_PI / area;
  const mRealType kalpha = 1.0 / (2.0*alpha);
  mRealType kmag, uk;

  Fk.resize(KList.kpts_cart.size());
  MaxKshell = KList.kshell.size() - 1;
  Fk_symm.resize(MaxKshell);

  for (int ks = 0, ki = 0; ks < Fk_symm.size(); ks++)
  {
    kmag = std::sqrt(KList.ksq[ki]);
    uk = knorm/kmag*(std::tanh(kmag*dgate)-std::erf(kalpha*kmag));
    Fk_symm[ks] = uk;
    while (ki < KList.kshell[ks + 1] && ki < Fk.size())
      Fk[ki++] = uk;
  }
}

EwaldScreen2D::mRealType EwaldScreen2D::evaluateLR_r0() const
{
  const mRealType kalpha = 1.0 / (2.0*alpha);
  mRealType kmax_d, kmax_a, kmax;
  // let 1-tanh(kd) < 5e-9
  kmax_d = 10.0/dgate;
  // let 1-erf(k/(2*alpha)) < 5e-9
  kmax_a = 8.3*alpha;
  kmax = std::max(kmax_d, kmax_a);
  // iFT[vlr_k] at r=0
  double dk = LR_kc/1e6;
  int nk = kmax/dk;
  double vlr_r0 = 0.0;
  for (int ik=0;ik<nk;ik++)
  {
    double k=ik*dk;
    vlr_r0 += std::tanh(k*dgate) - std::erf(k*kalpha);
  }
  vlr_r0 *= dk;
  return vlr_r0;
}

} // qmcplusplus
