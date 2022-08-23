/** @file GaussianOrbital.h
 *  * @brief Definition of member functions of Gaussian basis set
 *   */
#ifndef QMCPLUSPLUS_GAUSSIAN_ORBITAL
#define QMCPLUSPLUS_GAUSSIAN_ORBITAL

#include "QMCWaveFunctions/SPOSet.h"
#include "Particle/DistanceTableData.h"

namespace qmcplusplus
{

class GaussianOrbital : public SPOSet
{
public:
  GaussianOrbital(ParticleSet& target, ParticleSet& source, RealType cexpo, int ndim);
  ~GaussianOrbital();

  // phi[i][j] is phi_j(r_i), i.e. electron i in orbital j
  //  i \in [first, last)
  void evaluate_notranspose(
    const ParticleSet& P,
    int first,
    int last,
    ValueMatrix_t& phi,
    GradMatrix_t& dphi,
    ValueMatrix_t& d2phi) override;

  // plug r_i into all orbitals
  void evaluateVGL(
    const ParticleSet& P,
    int i,
    ValueVector_t& pvec,
    GradVector_t& dpvec,
    ValueVector_t& d2pvec
  ) override;
  void evaluateValue(const ParticleSet& P, int iat, ValueVector_t& pvec) override;

  // hessian matrix is needed by backflow
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& phi,
                            GradMatrix_t& dphi,
                            HessMatrix_t& d2phi_mat) override;

  // derivative of hessian is needed to optimize backflow
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& phi,
                            GradMatrix_t& dphi,
                            HessMatrix_t& d2phi_mat,
                            GGGMatrix_t& d3phi_mat) override;

  SPOSet* makeClone() const;
  // ---- begin required overrides
  void resetParameters(const opt_variables_type& optVariables) override {};
  void setOrbitalSetSize(int norbs) override {};
  // required overrides end ----
  void report(const std::string& pad) const override;
private:
  ParticleSet& targetPtcl;
  ParticleSet& sourcePtcl;
  RealType cexpo;
  const int ideitab;
  const int ndim;
  // helper functions
  const DistanceTableData::DistRow& getDistanceRow(const ParticleSet& P, const int i);
  const DistanceTableData::DisplRow& getDisplacementRow(const ParticleSet& P, const int i);
  //  phi(r=|{x, y, z}|)
  RealType operator()(const RealType rij);
  //  d/dx log( phi(r) ), etc.
  void gradient_log(GradType& dp, const RealType& rij, const PosType& drij);
  //  d^2/dx/dy log( phi(r) ), etc.
  void hessian_log(HessType& h, const GradType& dp);
  //  d^3/dx/dy/dx log( phi(r) ), etc.
  void gradHess_log(GGGType& g3, const GradType& dp);
  bool checkDerivatives;
};

} // qmcplusplus
#endif
