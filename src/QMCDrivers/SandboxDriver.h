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

#ifndef QMCPLUSPLUS_SANDBOXDRIVER_H
#define QMCPLUSPLUS_SANDBOXDRIVER_H
#include "QMCDrivers/QMCDriver.h"
#include "Particle/ParticleSetPool.h"
namespace qmcplusplus
{

class SandboxDriver : public QMCDriver
{
public:
  SandboxDriver(const ProjectData& project_data,
                     MCWalkerConfiguration& w,
                     TrialWaveFunction& psi,
                     QMCHamiltonian& h,
                     ParticleSetPool& ptclPool,
                     Communicate* comm);
  ~SandboxDriver() override final {};

  bool run() override final {return true;};
  bool put(xmlNodePtr q) override final {return true;}
  QMCRunType getRunType() override final { return QMCRunType::SANDBOX; }
private:
  ParticleSetPool& PtclPool;
  const size_t ndim;
};

} //qmcplusplus
#endif
