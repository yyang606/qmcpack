//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2017 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////



#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Numerics/OhmmsBlas.h"
#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"
#include "simd/simd.hpp"
#include "FakeSPO.h"


#include <stdio.h>
#include <string>

using std::string;


namespace qmcplusplus
{

typedef QMCTraits::ValueType ValueType;

template <typename T1, typename T2> void check_matrix(Matrix<T1> &a, Matrix<T2> &b)
{
  REQUIRE(a.size() == b.size());
  for (int i = 0; i < a.rows(); i++) {
    for (int j = 0; j < a.cols(); j++) {
      REQUIRE(a(i,j) == ValueApprox(b(i,j)));
    }
  }
}

TEST_CASE("DiracDeterminantBase_first", "[wavefunction][fermion]")
{
  FakeSPO *spo = new FakeSPO();
  spo->setOrbitalSetSize(3);
  DiracDeterminantBase ddb(spo);

  int norb = 3;
  ddb.set(0,norb);

  // occurs in call to registerData
  ddb.dpsiV.resize(norb);
  ddb.d2psiV.resize(norb);


  ParticleSet elec;

  elec.create(3);
  ddb.recompute(elec);

  Matrix<ValueType> b;
  b.resize(3,3);

  b(0,0) = 0.6159749342;
  b(0,1) = -0.2408954682;
  b(0,2) = -0.1646081192;
  b(1,0) = 0.07923894288;
  b(1,1) = 0.1496231042;
  b(1,2) = -0.1428117337;
  b(2,0) = -0.2974298429;
  b(2,1) = -0.04586322768;
  b(2,2) = 0.3927890292;

  check_matrix(ddb.psiM, b);


  DiracDeterminantBase::GradType grad;
  ValueType det_ratio = ddb.ratioGrad(elec, 0, grad);
  ValueType det_ratio1 = 0.178276269185;
  REQUIRE(det_ratio1 == ValueApprox(det_ratio));

  ddb.acceptMove(elec, 0);

  b(0,0) =  3.455170657;
  b(0,1) =  -1.35124809;
  b(0,2) = -0.9233316353;
  b(1,0) = 0.05476311768;
  b(1,1) = 0.1591951095;
  b(1,2) = -0.1362710138;
  b(2,0) = -2.235099338;
  b(2,1) = 0.7119205298;
  b(2,2) = 0.9105960265;

  check_matrix(ddb.psiM, b);


}

//#define DUMP_INFO

TEST_CASE("DiracDeterminantBase_second", "[wavefunction][fermion]")
{
  FakeSPO *spo = new FakeSPO();
  spo->setOrbitalSetSize(4);
  DiracDeterminantBase ddb(spo);

  int norb = 4;
  ddb.set(0,norb);

  // occurs in call to registerData
  ddb.dpsiV.resize(norb);
  ddb.d2psiV.resize(norb);


  ParticleSet elec;

  elec.create(4);
  ddb.recompute(elec);

  Matrix<ValueType> orig_a;
  orig_a.resize(4,4);
  orig_a = spo->a2;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < norb; j++) {
      orig_a(j,i) = spo->v2(i,j);
    }
  }

  //check_matrix(ddb.psiM, b);
  DiracMatrix<ValueType> dm;

  Matrix<ValueType> a_update1;
  a_update1.resize(4,4);
  a_update1 = spo->a2;
  for (int j = 0; j < norb; j++) {
    a_update1(j,0) = spo->v2(0,j);
  }

  Matrix<ValueType> a_update2;
  a_update2.resize(4,4);
  a_update2 = spo->a2;
  for (int j = 0; j < norb; j++) {
    a_update2(j,0) = spo->v2(0,j);
    a_update2(j,1) = spo->v2(1,j);
  }

  Matrix<ValueType> a_update3;
  a_update3.resize(4,4);
  a_update3 = spo->a2;
  for (int j = 0; j < norb; j++) {
    a_update3(j,0) = spo->v2(0,j);
    a_update3(j,1) = spo->v2(1,j);
    a_update3(j,2) = spo->v2(2,j);
  }


  DiracDeterminantBase::GradType grad;
  ValueType det_ratio = ddb.ratioGrad(elec, 0, grad);

  dm.invert(a_update1, true);
  ValueType det_update1 = dm.LogDet;
  ValueType log_ratio1 = det_update1 - ddb.LogValue;
  ValueType det_ratio1 = std::exp(log_ratio1);
#ifdef DUMP_INFO
  std::cout << "det 0 = " << std::exp(ddb.LogValue) << std::endl;
  std::cout << "det 1 = " << std::exp(det_update1) << std::endl;
  std::cout << "det ratio 1 = " << det_ratio1 << std::endl;
#endif
  //double det_ratio1 = 0.178276269185;

  REQUIRE(det_ratio1 == ValueApprox(det_ratio));

  ddb.acceptMove(elec, 0);


  ValueType det_ratio2 = ddb.ratioGrad(elec, 1, grad);
  dm.invert(a_update2, true);
  ValueType det_update2 = dm.LogDet;
  ValueType log_ratio2 = det_update2 - det_update1;
  ValueType det_ratio2_val = std::exp(log_ratio2);
#ifdef DUMP_INFO
  std::cout << "det 1 = " << std::exp(ddb.LogValue) << std::endl;
  std::cout << "det 2 = " << std::exp(det_update2) << std::endl;
  std::cout << "det ratio 2 = " << det_ratio2 << std::endl;
#endif
  //double det_ratio2_val = 0.178276269185;
  REQUIRE(std::abs(det_ratio2) == ValueApprox(det_ratio2_val));


  ddb.acceptMove(elec, 1);

  ValueType det_ratio3 = ddb.ratioGrad(elec, 2, grad);
  dm.invert(a_update3, true);
  ValueType det_update3 = dm.LogDet;
  ValueType log_ratio3 = det_update3 - det_update2;
  ValueType det_ratio3_val = std::exp(log_ratio3);
#ifdef DUMP_INFO
  std::cout << "det 2 = " << std::exp(ddb.LogValue) << std::endl;
  std::cout << "det 3 = " << std::exp(det_update3) << std::endl;
  std::cout << "det ratio 3 = " << det_ratio3 << std::endl;
#endif
  REQUIRE(det_ratio3 == ValueApprox(det_ratio3_val));
  //check_value(det_ratio3, det_ratio3_val);

  ddb.acceptMove(elec, 2);

  dm.invert(orig_a,false);

#ifdef DUMP_INFO
  std::cout << "original " << std::endl;
  std::cout << orig_a << std::endl;
  std::cout << "block update " << std::endl;
  std::cout << ddb.psiM << std::endl;
#endif

  check_matrix(orig_a, ddb.psiM);


}

}
