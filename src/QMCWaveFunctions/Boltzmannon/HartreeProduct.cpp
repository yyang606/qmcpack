#include <assert.h>
#include "QMCWaveFunctions/Boltzmannon/HartreeProduct.h"

namespace qmcplusplus
{

  HartreeProduct::HartreeProduct(SPOSetBasePtr const &spo, int first, int last) :
    Phi(spo), FirstIndex(first), LastIndex(last)
   ,SPOVGLTimer("HartreeProduct::spovgl",timer_level_fine)
  {
    OrbitalName = "HartreeProduct";
    registerTimers();

    // allocate matrices to hold SPO set data
    //  assume ground-state occupation i.e. number of orbitals = number of particles
    int nptcl = last-first;
    assert(nptcl>0);
    psiM.resize(nptcl,nptcl);
    dpsiM.resize(nptcl,nptcl);
    d2psiM.resize(nptcl,nptcl);
  }

  HartreeProduct::~HartreeProduct(){}
  void HartreeProduct::registerTimers()
  {
    TimerManager.addTimer(&SPOVGLTimer);
  }

  HartreeProduct::RealType
  HartreeProduct::evaluateLog(ParticleSet& P
    ,ParticleSet::ParticleGradient_t& G
    ,ParticleSet::ParticleLaplacian_t& L)
  {
    // YY: Is there an SPO function for just the diagonal of psiM etc.?
    SPOVGLTimer.start();
    Phi->evaluate(P,FirstIndex,LastIndex,psiM,dpsiM,d2psiM);
    SPOVGLTimer.stop();

    ValueType wf_val = 1.0;
    ValueType val;
    GradType  grad;
    ValueType lap;
    int nptcl = LastIndex-FirstIndex;
    for (int iptcl=0;iptcl<nptcl;iptcl++)
    {
      val  = psiM(iptcl,iptcl);
      wf_val *= val;
      grad = dpsiM(iptcl,iptcl)/val;
      lap  = d2psiM(iptcl,iptcl)/val-dot(grad,grad);
      G[FirstIndex+iptcl] += grad;
      L[FirstIndex+iptcl] += lap;
    }
    LogValue = evaluateLogAndPhase(wf_val,PhaseValue); // assign LogValue and PhaseValue
    return LogValue;
  }

}
