//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_RADIAL_JASTROW_BUILDER_H
#define QMCPLUSPLUS_RADIAL_JASTROW_BUILDER_H
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include <string>

namespace qmcplusplus
{

/** JastrowBuilder using an analytic 1d functor
 * Should be able to eventually handle all one and two body jastrows
 * although spline based ones will come later
 */


struct RadialJastrowBuilder: public WaveFunctionComponentBuilder
{
  RadialJastrowBuilder(ParticleSet& target, TrialWaveFunction& psi, PtclPoolType& psets);
  bool put(xmlNodePtr cur);
  PtclPoolType& ptclPool;

  template<typename RadfuncType>
  bool RadialJastrowBuilder::createJ2(xmlNodePtr cur, const std::string& jname)

};
}

#endif
