#include "WaveFunctionPlotter.h"

namespace qmcplusplus
{
WaveFunctionPlotter::WaveFunctionPlotter(xmlNodePtr cur,
                                         MCWalkerConfiguration& w,
                                         TrialWaveFunction& psi,
                                         QMCHamiltonian& h,
                                         Communicate* comm)
    : QMCDriver(w, psi, h, comm, "WaveFunctionPlotter")
{
}
bool WaveFunctionPlotter::run()
{
  app_log() << "Starting WF plotter\n";
  int nat = W.getTotalNum();
  app_log() << nat << "\n";
}
bool WaveFunctionPlotter::put(xmlNodePtr q)
{
}
} // namespace qmcplusplus
