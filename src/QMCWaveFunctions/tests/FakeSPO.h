#include "QMCWaveFunctions/SPOSetBase.h"

namespace qmcplusplus
{

typedef QMCTraits::ValueType ValueType;

class FakeSPO : public SPOSetBase
{
  public:

  Matrix<ValueType> a;
  Matrix<ValueType> a2;
  Vector<ValueType> v;
  Matrix<ValueType> v2;

  FakeSPO();
  virtual ~FakeSPO() {}

  virtual void report() {}
  virtual void resetParameters(const opt_variables_type& optVariables) {}
  virtual void resetTargetParticleSet(ParticleSet& P) {}
  virtual void setOrbitalSetSize(int norbs);

  virtual void
  evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

  virtual void
  evaluate(const ParticleSet& P, int iat,
           ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  virtual void
  evaluate(
    const ParticleSet& P,
    int iat,
    ValueVector_t& psi,
    GradVector_t& dpsi,
    HessVector_t& grad_grad_psi
  );

  virtual void evaluate_notranspose(
    const ParticleSet& P, int first, int last,
    ValueMatrix_t& logdet,
    GradMatrix_t& dlogdet,
    ValueMatrix_t& d2logdet
  );

};

}
