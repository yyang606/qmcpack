#include <QMCHamiltonians/ChargedSk.h>
#include <OhmmsData/AttributeSet.h>
#include <LongRange/KContainer.h>
#include <LongRange/StructFact.h>

namespace qmcplusplus
{

  ChargedSk::ChargedSk(ParticleSet& P)
      : Pinit(P)
  {
#ifndef USE_REAL_STRUCT_FACTOR
    APP_ABORT("ChargedSk: please recompile with USE_REAL_STRUCT_FACTOR=1");
#endif
    if(P.Lattice.SuperCellEnum==SUPERCELL_OPEN)
      APP_ABORT("ChargedSk is incompatible with open boundary conditions");

    // get particle information
    SpeciesSet& species = P.getSpeciesSet();
    int charge_idx = species.getAttribute("charge");
    nspecies = species.size();
    charge_vec.resize(nspecies);
    for(int s=0;s<nspecies;++s)
    { // record intra-species rho_k (rho_k^2 is automatically given in value_squared)
      charge_vec[s] = species(charge_idx,s);
    }
    ngroup = 2; // Sc(k), rhoc(k)
    reset();
  }

  void ChargedSk::reset()
  {
    myName   = "ChargedSk";
    UpdateMode.set(COLLECTABLE,1);
    ecut     = -1.0;
    nkpoints = -1;
  }

  
  QMCHamiltonianBase* ChargedSk::makeClone(ParticleSet& P, TrialWaveFunction& Psi)
  {
    return new ChargedSk(*this);
  }
 

  bool ChargedSk::put(xmlNodePtr cur)
  {
    using std::sqrt;
    reset();
    const k2_t& k2_init = Pinit.SK->KLists.ksq;

    std::string write_report = "no";
    OhmmsAttributeSet attrib;
    attrib.add(myName,"name");
    attrib.add(write_report,"report");
    attrib.put(cur);

    if(ecut<0.0)
      nkpoints = k2_init.size();
    else
    {
      RealType k2cut = 2.0*ecut;
      nkpoints = 0;
      for(int i=0;i<k2_init.size();++i)
        if(k2_init[i]<k2cut)
          nkpoints++;
    }
    crhok_r.resize(nkpoints);
    crhok_i.resize(nkpoints);
    
    if(nkpoints==0)
      APP_ABORT("ChargedSk::put  could not find any kpoints");
    
    ecut = .5*k2_init[nkpoints-1];

    if(write_report=="yes")
      report("  ");

    return true;
  }


  void ChargedSk::report(const std::string& pad)
  {
    app_log()<<pad<<"ChargedSk report"<< std::endl;
    app_log()<<pad<<"  name     = "<< myName   << std::endl;
    app_log()<<pad<<"  ecut     = "<< ecut     << std::endl;
    app_log()<<pad<<"  nkpoints = "<< nkpoints << std::endl;
    app_log()<<pad<<"end ChargedSk report"<< std::endl;
  }


  void ChargedSk::addObservables(PropertySetType& plist,BufferType& collectables)
  {
    myIndex=collectables.current();
    std::vector<RealType> tmp(2*2*nkpoints); // real & imag parts
    collectables.add(tmp.begin(),tmp.end());
  }


  void ChargedSk::registerCollectables(std::vector<observable_helper*>& h5desc, hid_t gid) const 
  {
    hid_t sgid=H5Gcreate(gid,myName.c_str(),0);
    observable_helper* oh=new observable_helper("kpoints");
    oh->open(sgid); // add to SkAll hdf group
    oh->addProperty(const_cast<std::vector<PosType>&>(Pinit.SK->KLists.kpts_cart),"value");
    h5desc.push_back(oh);
    std::vector<int> ng(2);
    ng[0] = 2;
    ng[1] = nkpoints;

    std::vector<std::string> group_name = {"csk","crhok"};
    for(int ig=0;ig<ngroup;++ig)
    {
      observable_helper* oh = new observable_helper(group_name[ig]);
      oh->set_dimensions(ng,myIndex+ig*2*nkpoints);
      oh->open(sgid);
      h5desc.push_back(oh);
    }
  }


  ChargedSk::Return_t ChargedSk::evaluate(ParticleSet& P)
  {
    // initialized
    std::fill(crhok_r.begin(),crhok_r.end(),0);
    std::fill(crhok_i.begin(),crhok_i.end(),0);
    RealType w=tWalker->Weight;

    // get references to density in reciprocal space
    const Matrix<RealType>& rhok_r = P.SK->rhok_r;
    const Matrix<RealType>& rhok_i = P.SK->rhok_i;

    // sum over species
    for (int k=0;k<nkpoints;k++)
    {
      for (int s=0;s<nspecies;s++)
      {
        crhok_r[k] += charge_vec[s]*rhok_r(s,k);
        crhok_i[k] += charge_vec[s]*rhok_i(s,k);
      }
    }

    // fill memory reserved in stat.h5
    int nkptot = rhok_r.cols();

    // record charged S(k)
    int ig     = 0;
    int kc = myIndex + ig*2*nkpoints;
    for(int k=0;k<nkpoints;++k,++kc)
      P.Collectables[kc] += w*(crhok_r[k]*crhok_r[k]);
    for(int k=0;k<nkpoints;++k,++kc)
      P.Collectables[kc] += w*(crhok_i[k]*crhok_i[k]);

    // record charged rho(k)
    ig = 1;
    kc = myIndex + ig*2*nkpoints;
    for(int k=0;k<nkpoints;++k,++kc)
      P.Collectables[kc] += w*crhok_r[k];
    for(int k=0;k<nkpoints;++k,++kc)
      P.Collectables[kc] += w*crhok_i[k];

    // guard memory integrity
    const int kc_max = myIndex + ngroup*2*nkpoints;
    if (kc > kc_max) APP_ABORT("ChargedSk is overwriting memory");

    return 0.0;
  }



}
