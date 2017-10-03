fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000   1) //////////////////////////////////////////////////////////////////////////////////////
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000   2) // This file is distributed under the University of Illinois/NCSA Open Source License.
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000   3) // See LICENSE file in top directory for details.
9a6ad0b7 (Jeongnim Kim   2006-08-08 17:02:07 +0000   4) //
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000   5) // Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000   6) //
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000   7) // File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
9f56aa0c (Anouar Benali  2016-09-27 23:45:29 +0000   8) //                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000   9) //                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
9f56aa0c (Anouar Benali  2016-09-27 23:45:29 +0000  10) //                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000  11) //                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
9f56aa0c (Anouar Benali  2016-09-27 23:45:29 +0000  12) //                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000  13) //                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000  14) //
9f56aa0c (Anouar Benali  2016-09-27 23:45:29 +0000  15) // File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000  16) //////////////////////////////////////////////////////////////////////////////////////
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000  17)     
fcaa478f (Anouar Benali  2016-09-27 00:41:58 +0000  18)     
9a6ad0b7 (Jeongnim Kim   2006-08-08 17:02:07 +0000  19) #include "QMCWaveFunctions/BasisSetFactory.h"
08f0d3d4 (Jeongnim Kim   2011-04-29 16:33:50 +0000  20) #include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  21) #include "QMCWaveFunctions/HarmonicOscillator/SHOSetBuilder.h"
f9758dd1 (Jeongnim Kim   2009-04-23 16:24:55 +0000  22) #if OHMMS_DIM == 3
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  23) #if !defined(QMC_COMPLEX)
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  24) #include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  25) #include "QMCWaveFunctions/MolecularOrbitals/GTOBuilder.h"
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  26) #include "QMCWaveFunctions/MolecularOrbitals/STOBuilder.h"
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  27) #include "QMCWaveFunctions/MolecularOrbitals/MolecularBasisBuilder.h"
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  28) #endif
f68e48ac (Jeongnim Kim   2011-07-20 15:46:17 +0000  29) 
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  30) #if defined(HAVE_EINSPLINE)
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  31) #include "QMCWaveFunctions/EinsplineSetBuilder.h"
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  32) #endif
c893984d (Ken Esler      2007-08-06 14:44:26 +0000  33) #endif
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000  34) #include "QMCWaveFunctions/CompositeSPOSet.h"
08f0d3d4 (Jeongnim Kim   2011-04-29 16:33:50 +0000  35) #include "QMCWaveFunctions/OptimizableSPOBuilder.h"
40155791 (Jeremy McMinis 2011-08-11 05:05:18 +0000  36) #include "QMCWaveFunctions/AFMSPOBuilder.h"
b4d43de9 (Jeongnim Kim   2008-06-09 19:00:32 +0000  37) #include "Utilities/ProgressReportEngine.h"
7cdc9ec5 (Jeongnim Kim   2008-09-19 17:49:20 +0000  38) #include "Utilities/IteratorUtility.h"
9a6ad0b7 (Jeongnim Kim   2006-08-08 17:02:07 +0000  39) #include "OhmmsData/AttributeSet.h"
9a6ad0b7 (Jeongnim Kim   2006-08-08 17:02:07 +0000  40) 
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  41) namespace qmcplusplus
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  42) {
9a6ad0b7 (Jeongnim Kim   2006-08-08 17:02:07 +0000  43) 
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  44)   //initialization of the static data of BasisSetFactory 
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  45)   std::map<std::string,BasisSetBuilder*> BasisSetFactory::basis_builders;
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  46)   BasisSetBuilder* BasisSetFactory::last_builder=0;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  47) 
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  48)   SPOSetBase* get_sposet(const std::string& name)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  49)   {
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  50)     int nfound = 0;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  51)     SPOSetBase* spo = 0;
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  52)     std::map<std::string,BasisSetBuilder*>::iterator it;
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  53)     for(it=BasisSetFactory::basis_builders.begin();
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  54)         it!=BasisSetFactory::basis_builders.end();++it)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  55)     {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  56)       std::vector<SPOSetBase*>& sposets = it->second->sposets;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  57)       for(int i=0;i<sposets.size();++i)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  58)       {
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  59)         SPOSetBase* sposet = sposets[i];
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  60)         if(sposet->objectName==name)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  61)         {
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  62)           spo = sposet;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  63)           nfound++;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  64)         }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  65)       }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  66)     }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  67)     if(nfound>1)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  68)     {
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  69)       write_basis_builders();
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  70)       APP_ABORT_TRACE(__FILE__, __LINE__, "get_sposet: requested sposet "+name+" is not unique");
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  71)     }
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  72)     //else if(spo==NULL)
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  73)     //{
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  74)     //  write_basis_builders();
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  75)     //  APP_ABORT("get_sposet: requested sposet "+name+" does not exist");
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  76)     //}
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  77)     return spo;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  78)   }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  79) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  80) 
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  81)   void write_basis_builders(const std::string& pad)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  82)   {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  83)     std::string pad2 = pad+"  ";
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  84)     std::map<std::string,BasisSetBuilder*>::iterator it;
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000  85)     for(it=BasisSetFactory::basis_builders.begin();it!=BasisSetFactory::basis_builders.end();++it)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  86)     {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  87)       const std::string& type = it->first;
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  88)       std::vector<SPOSetBase*>& sposets = it->second->sposets;
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  89)       app_log()<<pad<<"sposets for BasisSetBuilder of type "<<type<< std::endl;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  90)       for(int i=0;i<sposets.size();++i)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  91)       {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000  92)         app_log()<<pad2<<"sposet "<<sposets[i]->objectName<< std::endl;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  93)       }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  94)     }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  95)   }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  96) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  97) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000  98) 
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000  99) /** constructor
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 100)  * \param els reference to the electrons
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 101)  * \param psi reference to the wavefunction
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 102)  * \param ions reference to the ions
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 103)  */
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 104) BasisSetFactory::BasisSetFactory(ParticleSet& els, TrialWaveFunction& psi, PtclPoolType& psets):
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000 105)   OrbitalBuilderBase(els,psi), ptclPool(psets)
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 106) {
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 107)   ClassName="BasisSetFactory";
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 108) }
9a6ad0b7 (Jeongnim Kim   2006-08-08 17:02:07 +0000 109) 
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 110) BasisSetFactory::~BasisSetFactory()
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 111) {
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 112)   DEBUG_MEMORY("BasisSetFactory::~BasisSetFactory");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 113) }
9a6ad0b7 (Jeongnim Kim   2006-08-08 17:02:07 +0000 114) 
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 115) bool BasisSetFactory::put(xmlNodePtr cur)
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 116) {
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 117)   return true;
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 118) }
b4d43de9 (Jeongnim Kim   2008-06-09 19:00:32 +0000 119) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 120) BasisSetBuilder* BasisSetFactory::createBasisSet(xmlNodePtr cur,xmlNodePtr  rootNode)
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 121) {
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 122)   ReportEngine PRE(ClassName,"createBasisSet");
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 123)   std::string sourceOpt("ion0");
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 124)   std::string type("");
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 125)   std::string name("");
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 126)   std::string keyOpt("NMO"); //gaussian Molecular Orbital
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 127)   std::string transformOpt("yes"); //numerical Molecular Orbital
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 128)   std::string cuspC("no");  // cusp correction
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 129)   std::string cuspInfo("");  // file with precalculated cusp correction info
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 130)   OhmmsAttributeSet aAttrib;
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 131)   aAttrib.add(sourceOpt,"source");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 132)   aAttrib.add(cuspC,"cuspCorrection");
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 133)   aAttrib.add(type,"type");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 134)   aAttrib.add(keyOpt,"keyword");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 135)   aAttrib.add(keyOpt,"key");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 136)   aAttrib.add(name,"name");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 137)   aAttrib.add(transformOpt,"transform");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 138)   aAttrib.add(cuspInfo,"cuspInfo");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 139)   if(rootNode != NULL)
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 140)     aAttrib.put(rootNode);
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 141) 
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 142)   std::string type_in=type;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 143)   tolower(type);
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 144) 
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 145)   //when name is missing, type becomes the input
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 146)   if(name.empty()) name=type_in;
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 147) 
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 148)   BasisSetBuilder* bb=0;
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 149) 
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 150)   //check if builder can be reused
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 151)   std::map<std::string,BasisSetBuilder*>::iterator bbit=basis_builders.find(name);
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 152)   if(bbit!= basis_builders.end())
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 153)   {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 154)     app_log() << "Reuse BasisSetBuilder \""<<name << "\" type " << type_in << std::endl;
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 155)     app_log().flush();
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 156)     bb=(*bbit).second;
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 157)     bb->put(rootNode);
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 158)     return last_builder=bb;
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 159)   }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 160) 
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000 161)   //assign last_builder
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 162)   bb=last_builder;
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 163) 
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 164)   if (type == "composite")
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 165)   {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 166)     app_log() << "Composite SPO set with existing SPOSets." << std::endl;
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 167)     bb= new  CompositeSPOSetBuilder();
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 168)   }
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 169)   else if (type == "jellium" || type == "heg")
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 170)   {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 171)     app_log()<<"Electron gas SPO set"<< std::endl;
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 172)     bb = new ElectronGasBasisBuilder(targetPtcl,rootNode);
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 173)   }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 174)   else if (type == "sho")
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 175)   {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 176)     app_log()<<"Harmonic Oscillator SPO set"<< std::endl;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 177)     bb = new SHOSetBuilder(targetPtcl);
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 178)   }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 179)   else if (type == "linearopt")
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 180)   {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 181)     //app_log()<<"Optimizable SPO set"<< std::endl;
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 182)     bb = new OptimizableSPOBuilder(targetPtcl,ptclPool,rootNode);
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 183)   }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 184)   else if (type == "afm")
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 185)   {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 186)     //       app_log()<<"AFM SPO set"<< std::endl;
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 187)     bb = new AFMSPOBuilder(targetPtcl,ptclPool,rootNode);
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 188)   }
08f0d3d4 (Jeongnim Kim   2011-04-29 16:33:50 +0000 189) #if OHMMS_DIM ==3
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 190)   else if(type.find("spline")<type.size())
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 191)   {
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 192)     name=type_in;
f9758dd1 (Jeongnim Kim   2009-04-23 16:24:55 +0000 193) #if defined(HAVE_EINSPLINE)
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 194)     PRE << "EinsplineSetBuilder:  using libeinspline for B-spline orbitals.\n";
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 195)     bb = new EinsplineSetBuilder(targetPtcl,ptclPool,rootNode);
311d879a (Jeongnim Kim   2008-12-09 20:30:33 +0000 196) #else
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 197)     PRE.error("Einspline is missing for B-spline orbitals",true);
d86d2014 (Ken Esler      2007-08-06 17:23:08 +0000 198) #endif
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 199)   }
f68e48ac (Jeongnim Kim   2011-07-20 15:46:17 +0000 200) #if !defined(QMC_COMPLEX)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 201)   else if(type == "molecularorbital" || type == "mo")
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 202)   {
38112115 (Ye Luo         2017-07-11 20:10:51 -0500 203) #if defined(ENABLE_SOA)
ef98e419 (Ye Luo         2017-07-01 01:03:00 -0500 204)     PRE.error("Molecular orbital support is not ready on SoA builds. Stay tuned!",true);
ef98e419 (Ye Luo         2017-07-01 01:03:00 -0500 205) #endif
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 206)     ParticleSet* ions=0;
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 207)     //do not use box to check the boundary conditions
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 208)     if(targetPtcl.Lattice.SuperCellEnum==SUPERCELL_OPEN)
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 209)       targetPtcl.setBoundBox(false);
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 210)     //initialize with the source tag
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 211)     PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 212)     if(pit == ptclPool.end())
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 213)       PRE.error("Missing basisset/@source.",true);
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 214)     else
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 215)       ions=(*pit).second;
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 216)     if(transformOpt == "yes")
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 217)     {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 218)       app_log() << "Using MolecularBasisBuilder<NGOBuilder>" << std::endl;
4a2f1e15 (Miguel Morales 2010-04-19 19:26:09 +0000 219) #if QMC_BUILD_LEVEL>2
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 220)       bb = new MolecularBasisBuilder<NGOBuilder>(targetPtcl,*ions,cuspC=="yes",cuspInfo);
4a2f1e15 (Miguel Morales 2010-04-19 19:26:09 +0000 221) #else
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 222)       bb = new MolecularBasisBuilder<NGOBuilder>(targetPtcl,*ions,false);
4a2f1e15 (Miguel Morales 2010-04-19 19:26:09 +0000 223) #endif
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 224)     }
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 225)     else
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 226)     {
4a2f1e15 (Miguel Morales 2010-04-19 19:26:09 +0000 227) #if QMC_BUILD_LEVEL>2
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 228)       if(cuspC == "yes")
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 229)         app_log() <<" ****** Cusp Correction algorithm is only implemented in combination with numerical radial orbitals. Use transform=yes to enable this option. \n";
4a2f1e15 (Miguel Morales 2010-04-19 19:26:09 +0000 230) #endif
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 231)       if(keyOpt == "GTO")
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 232)         bb = new MolecularBasisBuilder<GTOBuilder>(targetPtcl,*ions);
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 233)       else if(keyOpt == "STO")
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 234)         bb = new MolecularBasisBuilder<STOBuilder>(targetPtcl,*ions);
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 235)     }
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 236)   }
f68e48ac (Jeongnim Kim   2011-07-20 15:46:17 +0000 237) #endif //!QMC_COMPLEX
f5921b47 (Paul Young     2017-05-17 09:47:49 -0500 238) #endif  //OHMMS_DIM==3
9c9057c2 (Paul Young     2017-05-16 11:18:32 -0500 239)   else
9c9057c2 (Paul Young     2017-05-16 11:18:32 -0500 240)   {
9c9057c2 (Paul Young     2017-05-16 11:18:32 -0500 241)     APP_ABORT("BasisSetFactory::createSPOSet cannot build basis set of unknown type "+type);
9c9057c2 (Paul Young     2017-05-16 11:18:32 -0500 242)   }
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 243)   PRE.flush();
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000 244) 
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000 245)   if(bb==0)
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000 246)     APP_ABORT_TRACE(__FILE__, __LINE__, "BasisSetFactory::createBasisSet\n  BasisSetBuilder creation failed.");
4f6b31db (Jeongnim Kim   2013-11-08 18:05:03 +0000 247) 
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 248)   if(bb == last_builder)
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 249)     app_log() << " Missing both \"@name\" and \"@type\". Use the last BasisSetBuilder." << std::endl;
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 250)   else
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 251)   {
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 252)     bb->setReportLevel(ReportLevel);
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 253)     bb->initCommunicator(myComm);
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 254)     bb->put(cur);
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 255)     app_log()<<"Built BasisSetBuilder \""<< name<< "\" of type "<< type << std::endl;
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 256)     basis_builders[name]=bb; //use name, if missing type is used
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 257)   }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 258)   last_builder = bb;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 259) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 260)   return bb;
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 261) }
9a6ad0b7 (Jeongnim Kim   2006-08-08 17:02:07 +0000 262) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 263) 
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 264) SPOSetBase* BasisSetFactory::createSPOSet(xmlNodePtr cur)
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 265) {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 266)   std::string bname("");
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 267)   std::string bsname("");
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 268)   std::string sname("");
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 269)   std::string type("");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 270)   OhmmsAttributeSet aAttrib;
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 271)   aAttrib.add(bname,"basisset");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 272)   aAttrib.add(bsname,"basis_sposet");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 273)   aAttrib.add(sname,"name");
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 274)   aAttrib.add(type,"type");
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 275)   //aAttrib.put(rcur);
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 276)   aAttrib.put(cur);
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 277) 
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 278)   //tolower(type);
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 279) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 280)   BasisSetBuilder* bb;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 281)   if(bname=="")
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 282)     bname=type;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 283)   if(type=="")
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 284)     bb = last_builder;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 285)   else if(basis_builders.find(type)!=basis_builders.end())
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 286)     bb = basis_builders[type];
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 287)   else
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 288)   {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 289)     std::string cname("");
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 290)     xmlNodePtr tcur=cur->children;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 291)     if(tcur!=NULL)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 292)       getNodeName(cname,cur);
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 293)     if(cname==basisset_tag)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 294)       bb = createBasisSet(tcur,cur);
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 295)     else
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 296)       bb = createBasisSet(cur,cur);
e5adf474 (Jeongnim Kim   2013-07-01 20:56:18 +0000 297)   }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 298)   if(bb)
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 299)   {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 300)     app_log()<<" Building SPOset "<<sname<<" with "<<bname<<" basis set."<< std::endl;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 301)     return bb->createSPOSet(cur);
9a6ad0b7 (Jeongnim Kim   2006-08-08 17:02:07 +0000 302)   }
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 303)   else
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 304)   {
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 305)     APP_ABORT("BasisSetFactory::createSPOSet Failed to create a SPOSet. basisBuilder is empty.");
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 306)     return 0;
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 307)   }
75c8ed1e (Jeremy McMinis 2013-04-26 00:14:53 +0000 308) }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 309) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 310) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 311) void BasisSetFactory::build_sposet_collection(xmlNodePtr cur)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 312) {
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 313)   xmlNodePtr parent = cur;
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 314)   std::string type("");
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 315)   OhmmsAttributeSet attrib;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 316)   attrib.add(type,"type");
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 317)   attrib.put(cur);
6ca6dda0 (Jeongnim Kim   2013-12-04 16:47:20 +0000 318)   //tolower(type); 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 319) 
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 320)   app_log()<<"building sposet collection of type "<<type<< std::endl;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 321) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 322)   BasisSetBuilder* bb = createBasisSet(cur,cur);
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 323)   xmlNodePtr element = parent->children;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 324)   int nsposets = 0;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 325)   while(element!=NULL)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 326)   {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 327)     std::string cname((const char*)(element->name));
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 328)     if(cname=="sposet")
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 329)     {
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 330)       std::string name("");
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 331)       OhmmsAttributeSet attrib;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 332)       attrib.add(name,"name");
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 333)       attrib.put(element);
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 334) 
fae28636 (Mark Berrill   2016-05-27 14:44:45 +0000 335)       app_log()<<"  Building SPOSet \""<<name<<"\" with "<<type<<" BasisSetBuilder"<< std::endl;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 336)       SPOSetBase* spo = bb->createSPOSet(element);
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 337)       spo->objectName = name;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 338)       nsposets++;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 339)     }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 340)     element = element->next;
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 341)   }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 342)   if(nsposets==0)
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 343)     APP_ABORT("BasisSetFactory::build_sposet_collection  no <sposet/> elements found");
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 344) }
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 345) 
adca8d2f (Jaron Krogel   2013-11-06 17:22:32 +0000 346) 
9a6ad0b7 (Jeongnim Kim   2006-08-08 17:02:07 +0000 347) }
