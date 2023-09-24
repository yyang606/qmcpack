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

#include "Particle/LongRange/Screen2D.h"

namespace qmcplusplus
{

Screen2D::Screen2D(ParticleSet& ref, mRealType kc_in)
  : LRHandlerBase(kc_in), dgate(ref.getLattice().dgate), mimg(ref.getLattice().mimg)
{
  if (ref.getLattice().ndim != 2)
    throw std::runtime_error("2D potential requires 2D Lattice");
  if (dgate < 0)
    throw std::runtime_error("screened potential requires distance_to_gate");
  LR_rc = ref.getLattice().LR_rc; // CoulombPBC needs get_rc() to createSpline4RbyVs
  LR_kc = ref.getLattice().LR_kc; // get_kc() is used in QMCFiniteSize
  area = ref.getLattice().Volume/ref.getLattice().R(2,2);
  // report
  app_log() << "    dgate = " << dgate << std::endl;
  app_log() << "    mimg  = " << mimg << std::endl;
  fillFk(ref.getSimulationCell().getKLists());
  llr = false;
}

void Screen2D::fillFk(const KContainer& KList)
{
  Fk.resize(KList.kpts_cart.size());
  MaxKshell = KList.kshell.size() - 1;
  Fk_symm.resize(MaxKshell);

  for (int ks = 0, ki = 0; ks < Fk_symm.size(); ks++)
  {
    Fk_symm[ks] = 0.0;
    while (ki < KList.kshell[ks + 1] && ki < Fk.size())
      Fk[ki++] = 0.0;
  }
}

Screen2D::mRealType Screen2D::evaluate(mRealType r, mRealType rinv) const
{
  mRealType vsr = 0.0;
  for (int m=-mimg;m<=mimg;m++)
  {
    double dm_p=2*m*2*dgate-2*dgate;
    double dm_n=2*m*2*dgate;
    vsr += -1/std::sqrt(r*r + dm_p*dm_p) + 1/std::sqrt(r*r + dm_n*dm_n);
  }
  return vsr;
}

} // qmcplusplus
