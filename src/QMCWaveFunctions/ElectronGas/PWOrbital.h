/** @file PWOrbital.h
 *  * @brief Definition of member functions of plane wave (PW) basis set
 *   */
#ifndef QMCPLUSPLUS_PW_ORBITAL
#define QMCPLUSPLUS_PW_ORBITAL

#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{

class PWOrbital : public SPOSet
{
public:
  PWOrbital(const std::vector<PosType>& kpts_cart);
  ~PWOrbital();

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

  // ---- begin required overrides
  SPOSet* makeClone() const {return new PWOrbital(*this);}
  void resetParameters(const opt_variables_type& optVariables) override {APP_ABORT("not implemented")};
  void setOrbitalSetSize(int norbs) override {APP_ABORT("not implemented")};
  // required overrides end ----
  void report(const std::string& pad) const override;
private:
  std::vector<PosType> K; // K vectors
  std::vector<RealType> mK2; // minus K^2
};

} // qmcplusplus
#endif
