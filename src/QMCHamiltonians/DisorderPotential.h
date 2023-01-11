#ifndef QMCPLUSPLUS_DISORDER_POTENTIAL_H
#define QMCPLUSPLUS_DISORDER_POTENTIAL_H

//#include <algorithm>
#include "QMCHamiltonians/OperatorBase.h"

namespace qmcplusplus
{
class DisorderPotential : public OperatorBase
{
public:
  DisorderPotential(const ParticleSet& ions, ParticleSet& P)
   : itab(P.addTable(ions)), nsite(ions.getTotalNum())
  {
    setEnergyDomain(POTENTIAL);
    oneBodyQuantumDomain(P);
    widths.resize(nsite);
    heights.resize(nsite);
    std::fill(widths.begin(), widths.end(), 5.0);
    std::fill(heights.begin(), heights.end(), 1.0);
  };
  ~DisorderPotential(){};
  Return_t evaluate(ParticleSet& P) override;

  bool put(xmlNodePtr cur) override;
  bool get(std::ostream& os) const override;
  std::unique_ptr<OperatorBase> makeClone(ParticleSet& P, TrialWaveFunction& psi) override;
  // ---- begin required overrides
  void resetTargetParticleSet(ParticleSet& P) override {APP_ABORT("not implemented");};
  std::string getClassName() const override {return "disorder";};
  // required overrides end ----
private:
  std::vector<RealType> widths, heights;
  const int itab;
  const size_t nsite;
};
} // qmcplusplus
#endif
