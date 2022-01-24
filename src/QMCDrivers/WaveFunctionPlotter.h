#ifndef QMCPLUSPLUS_WAVEFUNCTIONPLOT_H
#define QMCPLUSPLUS_WAVEFUNCTIONPLOT_H

#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

/** Plot TrialWaveFunction values, gradients and laplacians
*/
class WaveFunctionPlotter : public QMCDriver
{
public:
  /// type definition
  using LogValueType = WaveFunctionComponent::LogValueType;

  /// Constructor.
  WaveFunctionPlotter(xmlNodePtr cur,
                      MCWalkerConfiguration& w,
                      TrialWaveFunction& psi,
                      QMCHamiltonian& h,
                      Communicate* comm);

  ~WaveFunctionPlotter(){};

  bool run();
  bool put(xmlNodePtr q);
  QMCRunType getRunType() { return QMCRunType::WF_PLOT; }
private:
  int nstep, ipart, jpart;
};

} // namespace qmcplusplus
#endif
