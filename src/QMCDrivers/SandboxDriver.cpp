#include "QMCDrivers/SandboxDriver.h"

namespace qmcplusplus
{

SandboxDriver::SandboxDriver(const ProjectData& project_data,
                                   MCWalkerConfiguration& w,
                                   TrialWaveFunction& psi,
                                   QMCHamiltonian& h,
                                   ParticleSetPool& ptclPool,
                                   Communicate* comm)
  : QMCDriver(project_data, w, psi, h, comm, "SandboxDriver"),
  PtclPool(ptclPool),
  ndim(w.getLattice().ndim)
{
}


}
