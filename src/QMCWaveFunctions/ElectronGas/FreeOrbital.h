//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_FREE_ORBITAL
#define QMCPLUSPLUS_FREE_ORBITAL

#include "QMCWaveFunctions/SPOSet.h"

namespace qmcplusplus
{
class FreeOrbital : public SPOSet
{
public:
  FreeOrbital(const std::string& my_name, const std::vector<PosType>& kpts_cart);
  ~FreeOrbital();

  // phi[i][j] is phi_j(r_i), i.e. electron i in orbital j
  //  i \in [first, last)
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix_t& phi,
                            GradMatrix_t& dphi,
                            ValueMatrix_t& d2phi) override;

  // plug r_i into all orbitals
  void evaluateVGL(const ParticleSet& P, int i, ValueVector_t& pvec, GradVector_t& dpvec, ValueVector_t& d2pvec) override;
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

  void report(const std::string& pad) const override;
  // ---- begin required overrides
  void resetParameters(const opt_variables_type& optVariables) override {};
  SPOSet* makeClone() const override { return new FreeOrbital(*this); }
  void setOrbitalSetSize(int norbs) override { throw std::runtime_error("not implemented"); }
  // required overrides end ----
private:
  const std::vector<PosType> kvecs; // kvecs vectors
  const int mink;                   // minimum k index
  const int maxk;                   // maximum number of kvecs vectors
  std::vector<RealType> k2neg;      // minus kvecs^2
};

} // namespace qmcplusplus
#endif
