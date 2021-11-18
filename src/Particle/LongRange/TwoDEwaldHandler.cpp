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


#include "TwoDEwaldHandler.h"

namespace qmcplusplus
{
void TwoDEwaldHandler::initBreakup(ParticleSet& ref)
{
  LR_rc = ref.Lattice.LR_rc;
  LR_kc = ref.Lattice.LR_kc;
  Sigma = ref.Lattice.alpha;
  app_log() << "  Sigma=" << Sigma << std::endl;
  Volume = ref.Lattice.Volume/ref.Lattice.R(2, 2);
  app_log() << "  Area=" << Volume << std::endl;
  fillFk(ref.SK->KLists);
}

TwoDEwaldHandler::TwoDEwaldHandler(const TwoDEwaldHandler& aLR, ParticleSet& ref)
    : LRHandlerBase(aLR), Sigma(aLR.Sigma), Volume(aLR.Volume)
{}

void TwoDEwaldHandler::fillFk(KContainer& KList)
{
  Fk.resize(KList.kpts_cart.size());
  const std::vector<int>& kshell(KList.kshell);
  MaxKshell = kshell.size() - 1;
  Fk_symm.resize(MaxKshell);
  kMag.resize(MaxKshell);
  mRealType knorm           = 2.0 * M_PI / Volume;
  mRealType oneovertwosigma = 1.0 / (2.0 * Sigma);
  int nxy, nz;
  for (int ks = 0, ki = 0; ks < Fk_symm.size(); ks++)
  {
    kMag[ks]    = std::sqrt(KList.ksq[ki]);
    mRealType uk = knorm * erfc(kMag[ks] * oneovertwosigma) / kMag[ks];
    nxy = nz = 0;
    while (ki < KList.kshell[ks + 1] && ki < Fk.size())
    {
      if (std::abs(KList.kpts[ki][2]) < 1e-8)
      {
        Fk[ki] = uk;
        nxy++;
      } else {
        Fk[ki] = 0;
        nz++;
      }
      //app_log() << KList.kpts_cart[ki] << " " << kMag[ks] << " " << Fk[ki] << std::endl;
      ki++;
    }
    Fk_symm[ks] = uk*nxy/(nxy+nz);
    //app_log() << kMag[ks] << " " << uk << " " << nxy << " " << nz << " " << Fk_symm[ks]  << std::endl;
    //       app_log()<<kMag[ks]<<" "<<uk<< std::endl;
  }
  app_log().flush();
}
} // namespace qmcplusplus
