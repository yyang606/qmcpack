//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//
// File created by: Yubo "Paul" Yang, yubo.paul.yang@gmail.com, CCQ @ Flatiron
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_GAUSSIAN_ORBITAL
#define QMCPLUSPLUS_GAUSSIAN_ORBITAL

#include "QMCWaveFunctions/SPOSet.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

class GaussianOrbital : public SPOSet
{
public:
  GaussianOrbital(const std::string& my_name, ParticleSet& target, ParticleSet& source, RealType cexpo, const TinyVector<int, OHMMS_DIM>& pbc_images);
  ~GaussianOrbital();

  // phi[i][j] is phi_j(r_i), i.e. electron i in orbital j
  //  i \in [first, last)
  void evaluate_notranspose(
    const ParticleSet& P,
    int first,
    int last,
    ValueMatrix& phi,
    GradMatrix& dphi,
    ValueMatrix& d2phi) override;

  // plug r_i into all orbitals
  void evaluateVGL(
    const ParticleSet& P,
    int i,
    ValueVector& pvec,
    GradVector& dpvec,
    ValueVector& d2pvec
  ) override;
  void evaluateValue(const ParticleSet& P, int iat, ValueVector& pvec) override;

  // hessian matrix is needed by backflow
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& phi,
                            GradMatrix& dphi,
                            HessMatrix& d2phi_mat) override;

  // derivative of hessian is needed to optimize backflow
  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& phi,
                            GradMatrix& dphi,
                            HessMatrix& d2phi_mat,
                            GGGMatrix& d3phi_mat) override;

  void report(const std::string& pad) const override;
  std::string getClassName() const override { return "GaussianOrbital"; }
  // ---- begin required overrides
  std::unique_ptr<SPOSet> makeClone() const override { return std::make_unique<GaussianOrbital>(*this); }
  void setOrbitalSetSize(int norbs) override {};
  // required overrides end ----
private:
  ParticleSet& targetPtcl;
  ParticleSet& sourcePtcl;
  RealType cexpo;
  const int ideitab;
  const int ndim;
  // helper functions
  const DistanceTable::DistRow& getDistanceRow(const ParticleSet& P, const int i);
  const DistanceTable::DisplRow& getDisplacementRow(const ParticleSet& P, const int i);
  //  phi(r=|{x, y, z}|)
  RealType operator()(const RealType rij);
  //  d/dx log( phi(r) ), etc.
  void gradient_log(GradType& dp, const RealType& rij, const PosType& drij);
  //  d^2/dx/dy log( phi(r) ), etc.
  void hessian_log(HessType& h, const GradType& dp);
  //  d^3/dx/dy/dx log( phi(r) ), etc.
  void gradHess_log(GGGType& g3, const GradType& dp);
  bool checkDerivatives;
  // Number of Cell images for the evaluation of the orbital with PBC. If No PBC, should be 0;
  const TinyVector<int, OHMMS_DIM> PBCImages;
  inline int indexPBCImage(const int i) const {return ((i%2)*2-1)*((i+1)/2);};
};

} // qmcplusplus
#endif
