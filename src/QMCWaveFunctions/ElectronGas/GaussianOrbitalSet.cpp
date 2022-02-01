#include "GaussianOrbitalSet.h"

namespace qmcplusplus
{

GaussianOrbitalSet::GaussianOrbitalSet(RealType cexpo_in)
 : cexpo(cexpo_in) {}

GaussianOrbitalSet::~GaussianOrbitalSet(){}

void GaussianOrbitalSet::evaluate_notranspose(
  const ParticleSet& P,
  int first,
  int last,
  ValueMatrix_t& logdet,
  GradMatrix_t& dlogdet,
  ValueMatrix_t& d2logdet)
{
}

void GaussianOrbitalSet::report(const std::string& pad) const
{
  app_log() << pad << "GaussianOrbital report" << std::endl;
  app_log() << pad << "  cexpo    = " << cexpo << std::endl;
  app_log() << pad << "end GaussianOrbital report" << std::endl;
  app_log().flush();
}
} // qmcplusplus
