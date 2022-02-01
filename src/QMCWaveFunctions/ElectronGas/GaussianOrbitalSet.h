/** @file GaussianOrbitalSet.h
 *  * @brief Definition of member functions of Gaussian basis set
 *   */
#ifndef QMCPLUSPLUS_GAUSSIAN_ORBITAL_SET
#define QMCPLUSPLUS_GAUSSIAN_ORBITAL_SET

#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{

class GaussianOrbitalSet : public SPOSet
{
public:
  GaussianOrbitalSet();
  ~GaussianOrbitalSet();
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& logdet,
                            GradMatrix_t& dlogdet,
                            ValueMatrix_t& d2logdet) override;
  // ---- begin required overrides
  void resetParameters(const opt_variables_type& optVariables) override {};
  void setOrbitalSetSize(int norbs) override {};
  void evaluateValue(const ParticleSet& P, int iat, ValueVector_t& psi) override {};
  void evaluateVGL(
    const ParticleSet& P,
    int iat,
    ValueVector_t& psi,
    GradVector_t& dpsi,
    ValueVector_t& d2psi
  ) override {};
  // required overrides end ----
};

} // qmcplusplus
#endif
