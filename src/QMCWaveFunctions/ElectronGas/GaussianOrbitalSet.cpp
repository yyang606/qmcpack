#include "GaussianOrbitalSet.h"

namespace qmcplusplus
{

GaussianOrbitalSet::GaussianOrbitalSet(){}
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

} // qmcplusplus
