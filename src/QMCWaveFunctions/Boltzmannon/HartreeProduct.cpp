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
    int nptcl = last-first;
    assert(nptcl>0);
    int norb = nptcl; //  assume ground-state occupation i.e. number of orbitals = number of particles
    psiM.resize(nptcl,norb);
    dpsiM.resize(nptcl,norb);
    d2psiM.resize(nptcl,norb);
  }

  OrbitalBasePtr HartreeProduct::makeClone(ParticleSet& qp) const
  { // need to clone SPO set and point it to new particle set qp
    SPOSetBasePtr spo = Phi->makeClone();
    spo->resetTargetParticleSet(qp);
    HartreeProduct* myclone = new HartreeProduct(spo,FirstIndex,LastIndex);
    return myclone;
  }
  
  void HartreeProduct::resetTargetParticleSet(ParticleSet& qp)
  {
    Phi->resetTargetParticleSet(qp);
  }

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
      grad = dpsiM(iptcl,iptcl)/val; // grad_psi_over_psi
      lap  = d2psiM(iptcl,iptcl)/val-dot(grad,grad); // lap_psi_over_psi

      // accumulate VGL changes
      wf_val *= val;
      G[FirstIndex+iptcl] += grad;
      L[FirstIndex+iptcl] += lap;
    }
    LogValue = evaluateLogAndPhase(wf_val,PhaseValue); // assign LogValue and PhaseValue
    return LogValue;
  }

}
