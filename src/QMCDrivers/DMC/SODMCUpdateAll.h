#ifndef QMCPLUSPLUS_SODMC_ALLPARTICLE_UPDATE_H
#define QMCPLUSPLUS_SODMC_ALLPARTICLE_UPDATE_H
#include "QMCDrivers/QMCUpdateBase.h"
namespace qmcplusplus
{
class SODMCUpdateAll : public QMCUpdateBase
{
public:
  /// Constructor.
  SODMCUpdateAll(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, RandomGenerator& rg);
  ///destructor
  ~SODMCUpdateAll() override;

  void advanceWalker(Walker_t& thisWalker, bool recompute) override;

  /// Copy Constructor (disabled)
  SODMCUpdateAll(const SODMCUpdateAll&) = delete;
  /// Copy operator (disabled).
  SODMCUpdateAll& operator=(const SODMCUpdateAll&) = delete;
};
} // namespace qmcplusplus

#endif
