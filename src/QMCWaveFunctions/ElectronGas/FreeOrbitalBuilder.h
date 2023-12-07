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

#ifndef QMCPLUSPLUS_FREE_ORBITAL_BUILDER_H
#define QMCPLUSPLUS_FREE_ORBITAL_BUILDER_H

#include "QMCWaveFunctions/SPOSetBuilder.h"

namespace qmcplusplus
{
class FreeOrbitalBuilder : public SPOSetBuilder
{
public:
  FreeOrbitalBuilder(ParticleSet& els, Communicate* comm, xmlNodePtr cur);
  ~FreeOrbitalBuilder() {}

  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;
  double calc_kf(const PtclOnLatticeTraits::ParticleLayout lattice, const size_t nelec);

private:
  ParticleSet& targetPtcl;
  bool in_list(const int j, const std::vector<int> l);
};
} // namespace qmcplusplus
#endif
