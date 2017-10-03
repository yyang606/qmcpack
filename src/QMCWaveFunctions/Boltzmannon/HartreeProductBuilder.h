#ifndef QMCPLUSPLUS_HARTREE_PRODUCT_BUILDER_H
#define QMCPLUSPLUS_HARTREE_PRODUCT_BUILDER_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/BasisSetFactory.h"

namespace qmcplusplus
{

class HartreeProductBuilder: public OrbitalBuilderBase
{

  public:
    HartreeProductBuilder(ParticleSet& p, TrialWaveFunction& psi);
    bool put(xmlNodePtr cur); // parse xml node pointed to by cur (should be <hartree_product>)

};

}

#endif // QMCPLUSPLUS_HARTREE_PRODUCT_BUILDER_H
