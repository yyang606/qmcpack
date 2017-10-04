#ifndef QMCPLUSPLUS_HARTREEPRODUCT_H
#define QMCPLUSPLUS_HARTREEPRODUCT_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Utilities/NewTimer.h"

namespace qmcplusplus
{

/** wavefunction value = product of single-particle orbitals (SPO)
 * 
 * The methods in this class construct a many-body wavefunction for N particles by taking the
 *  product of N SPOs.
 *
 * The main data structure stored in this class is 'Phi', which is a pointer to an SPOSetBase. 
 *  When Phi is evaluated at a particle configuration, results are stored in psiM,dpsiM,d2psiM.
 *
 */
class HartreeProduct: public OrbitalBase
{
public:
  /* Constructor
   *  make space in psiM,dpsiM,d2psiM to hold SPO data
   *  also initialize variables and timers
   *
   *  spo: the set of SPOs to make a product.
   *  first: index of the first particle to go into the SPOs in 
   *   global quantum (a.k.a. target) particle set.
   *  last: index of the last particle.
   */
  HartreeProduct(SPOSetBasePtr const &spo, int first, int last);

  /* Destructor
   *  do nothing
   *  assume psiM, dpsiM, d2psiM classes have proper destructors
   *  assume SPO set is deallocated elsewhere 
   *   (because SPO set is not allocated in this class)
   */
  ~HartreeProduct(){};

  /* Clone method for OpenMP
   *  SPO set needs to be cloned?
   */
  OrbitalBasePtr makeClone(ParticleSet& qp) const;

  /* evaluateLog method
   * evaluate wavefunction value, (particle) gradient, (particle) laplacian
   *  this is the most basic functionality of a wavefunction. This function:
   *   1. calculates and stores LogValue and PhaseValue
   *   2. calculates grdient and store in G
   *   3. calculates laplaciant and store in L
   */
  RealType
  evaluateLog(ParticleSet& P,
              ParticleSet::ParticleGradient_t& G,
              ParticleSet::ParticleLaplacian_t& L);

  /* evaluate method
   *  same as evaluateLog but returns wavefunction value instead of its log
   */
  ValueType 
  evaluate(ParticleSet& P,
           ParticleSet::ParticleGradient_t& G,
           ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }

  // SPO set 'Phi' is hard-wired to a particle set at build time,
  //  thus Phi needs to be updated if particle set changes
  void resetTargetParticleSet(ParticleSet& P); 
  void registerTimers();

  // ---- override pure virtual functions using vacuous functions ---- //
  // no wavefunction parameter to optimize
  void checkInVariables(opt_variables_type& o){};
  void checkOutVariables(const opt_variables_type& o){};
  void resetParameters(const opt_variables_type& active){};
  void reportStatus(std::ostream& os){};
  // no support for particle-by-particle move
  ValueType ratio(ParticleSet& P,int iat)
  {
   APP_ABORT("ratio(P,iat) not implemented");
  }
  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL)
  {
   APP_ABORT("ratio(P,dG,dL) not implemented");
  }
  RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
  {
   APP_ABORT("evaluateLog(P,bug) not implemented");
  }
  RealType registerData(ParticleSet& P, PooledData<RealType>& buf)
  {
   APP_ABORT("registerData not implemented");
  }
  RealType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false)
  {
   APP_ABORT("updateBuffer not implemented");
  }
  void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
   APP_ABORT("copyFromBuffer not implemented");
  }
  void acceptMove(ParticleSet& P, int iat)
  {
   APP_ABORT("acceptMove not implemented");
  }
  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat)
  {
   APP_ABORT("update(P,dG,dL,iat) not implemented");
  }
  void restore(int iat)
  {
   APP_ABORT("restore not implemented");
  }

  // ---- accessors and private variables ---- //
  SPOSetBasePtr getPhi()
  {
    return Phi;
  }

private:
  int FirstIndex, LastIndex; // index of first and last particles in global quantum particle set
  SPOSetBasePtr Phi;         // pointer to single-particle orbital (SPO) set
  ValueMatrix_t psiM;        // values of SPOs (psiM[i,j]=Phi_j(r_i))
  GradMatrix_t  dpsiM;       // particle gradients of SPOs
  ValueMatrix_t d2psiM;      // particle laplacians of SPOs
  NewTimer SPOVGLTimer;      // timer for evaluating the value, gradient, and laplacian (VGL) of SPOs
};

}
#endif //QMCPLUSPLUS_HARTREEPRODUCT_H
