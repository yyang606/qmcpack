#ifndef QMCPLUSPLUS_WAVEFUNCTIONPLOT_H
#define QMCPLUSPLUS_WAVEFUNCTIONPLOT_H

#include "QMCDrivers/QMCDriver.h"
namespace qmcplusplus
{

/** Plot TrialWaveFunction values, gradients and laplacians
 * or external potential
*/
class WaveFunctionPlotter : public QMCDriver
{
public:
  /// type definition
  using LogValueType = WaveFunctionComponent::LogValueType;

  /// default overrides
  WaveFunctionPlotter(const ProjectData& project_data,
                      MCWalkerConfiguration& w,
                      TrialWaveFunction& psi,
                      QMCHamiltonian& h,
                      Communicate* comm);

  ~WaveFunctionPlotter(){};
  bool put(xmlNodePtr q) override;
  QMCRunType getRunType() override { return QMCRunType::WF_PLOT; }

  bool run() override;
private:
  std::string extpot;
  /// functionalities
  void plotExternalPotential();
};

} // namespace qmcplusplus
#endif
