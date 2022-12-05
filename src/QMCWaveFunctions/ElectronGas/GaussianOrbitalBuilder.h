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

#ifndef QMCPLUSPLUS_GAUSSIAN_ORBITAL_BUILDER_H
#define QMCPLUSPLUS_GAUSSIAN_ORBITAL_BUILDER_H

#include "QMCWaveFunctions/SPOSetBuilder.h"

namespace qmcplusplus
{
class GaussianOrbitalBuilder : public SPOSetBuilder
{
public:
  using PSetMap = std::map<std::string, const std::unique_ptr<ParticleSet>>;
  GaussianOrbitalBuilder(ParticleSet& els, const PSetMap& psets, Communicate* comm, xmlNodePtr cur);
  ~GaussianOrbitalBuilder();

  std::unique_ptr<SPOSet> createSPOSetFromXML(xmlNodePtr cur) override;
private:
  ParticleSet& targetPtcl;
  std::string sourceName;
  ///reference to the particleset pool
  const PSetMap& particleSets;
};
} // qmcplusplus
#endif
