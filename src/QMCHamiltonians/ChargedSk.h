#ifndef QMCPLUSPLUS_CHARGED_SK_H
#define QMCPLUSPLUS_CHARGED_SK_H

#include <QMCHamiltonians/QMCHamiltonianBase.h>

namespace qmcplusplus
{

class ChargedSk : public QMCHamiltonianBase
{
 public:
 
  typedef std::vector<RealType> k2_t;
  typedef std::vector<PosType>  pts_t;

  //data members
  int nspecies, ngroup;
  RealType ecut;
  int nkpoints;
  const ParticleSet& Pinit;
  std::vector<RealType> charge_vec; // charge of species
  std::vector<RealType> crhok_r, crhok_i; // species-summed charged rho(k)

  //constructor/destructor
  ChargedSk(ParticleSet& P);
  ~ChargedSk() { }

  //standard interface
  QMCHamiltonianBase* makeClone(ParticleSet& P, TrialWaveFunction& psi);
  bool put(xmlNodePtr cur);
  Return_t evaluate(ParticleSet& P);
  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P); 
  }

  //required for Collectables interface
  void addObservables(PropertySetType& plist,BufferType& olist);
  void registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const ;

  //should be empty for Collectables interface
  void resetTargetParticleSet(ParticleSet& P)                      { }
  void setObservables(PropertySetType& plist)                      { }
  void setParticlePropertyList(PropertySetType& plist, int offset) { }

#if !defined(REMOVE_TRACEMANAGER)
  void checkout_scalar_arrays(TraceManager& tm)                    { }
  void collect_scalar_samples()                                    { }
  void delete_scalar_arrays()                                      { }
#endif

  //obsolete?
  bool get(std::ostream& os) const { return false; }

  //local functions
  void reset();
  void report(const std::string& pad);

};

}

#endif
