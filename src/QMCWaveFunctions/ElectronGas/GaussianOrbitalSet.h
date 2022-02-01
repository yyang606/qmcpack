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
  GaussianOrbitalSet(ParticleSet& target, ParticleSet& source, RealType cexpo, int ndim);
  ~GaussianOrbitalSet();

  // phi[i][j] is phi_j(r_i), i.e. electron i in orbital j
  //  i \in [first, last]
  void evaluate_notranspose(
    const ParticleSet& P,
    int first,
    int last,
    ValueMatrix_t& phi,
    GradMatrix_t& dphi,
    ValueMatrix_t& d2phi) override;

  // ---- begin required overrides
  void resetParameters(const opt_variables_type& optVariables) override {};
  void setOrbitalSetSize(int norbs) override {};
  void evaluateValue(const ParticleSet& P, int iat, ValueVector_t& psi) override {};
  void evaluateVGL(
    const ParticleSet& P,
    int i,
    ValueVector_t& pvec,
    GradVector_t& dpvec,
    ValueVector_t& d2pvec
  ) override {};
  // required overrides end ----
  void report(const std::string& pad) const override;
private:
  ParticleSet& targetPtcl;
  ParticleSet& sourcePtcl;
  RealType cexpo;
  const int ideitab;
  const int ndim;
};

} // qmcplusplus
#endif
