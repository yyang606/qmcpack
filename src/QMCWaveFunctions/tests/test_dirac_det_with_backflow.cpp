#include "Message/catch_mpi_main.hpp"
//#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "FakeSPO.h"

#include "QMCWaveFunctions/Fermion/DiracDeterminantWithBackflow.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Fermion/Backflow_ee.h"

using std::string;
namespace qmcplusplus
{

typedef QMCTraits::RealType   RealType;
typedef QMCTraits::ValueType ValueType;

TEST_CASE("DiracDeterminantWithBackflow_nobf", "[wavefunction][fermion]")
{

  ParticleSet elec;
  elec.setName("e");
  elec.create(3);

  double rcut = 9.90104015;

  // make a spline that evaluates to 0 constant
  BsplineFunctor<RealType> bsp;
  bsp.resize(1);
  bsp.cutoff_radius = rcut;
  
  // make a e-e backflow transformation that does nothing
  Backflow_ee<BsplineFunctor<RealType>> tbf(elec, elec);
  tbf.addFunc(0, 0, &bsp);

  // make backflow transform class
  BackflowTransformation bftrans(elec);
  bftrans.cutOff = rcut;
  bftrans.bfFuns.push_back((BackflowFunctionBase*)&tbf);

  // make determinant
  FakeSPO *spo = new FakeSPO();
  spo->setOrbitalSetSize(3);
  DiracDeterminantWithBackflow det(elec, spo, &bftrans);

}

}
