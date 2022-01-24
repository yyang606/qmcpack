#include "WaveFunctionPlotter.h"

namespace qmcplusplus
{
WaveFunctionPlotter::WaveFunctionPlotter(xmlNodePtr cur,
                                         MCWalkerConfiguration& w,
                                         TrialWaveFunction& psi,
                                         QMCHamiltonian& h,
                                         Communicate* comm)
    : QMCDriver(w, psi, h, comm, "WaveFunctionPlotter")
{
}
bool WaveFunctionPlotter::run()
{
  app_log() << "Starting WF plotter\n";
  int npart = W.getTotalNum();
  RealType dr = 0.01, r, eloc;
  ValueType logpsi;
  if (jpart < 0) jpart = npart+jpart;
  std::string fname = "test.dat";
  std::ofstream eout(fname.c_str());
  app_log() << "Moving " << jpart << " @ " << W.R[jpart] << "\n"
            << " onto "  << ipart << " @ " << W.R[ipart] << " " << npart << "\n";

  W.loadWalker(**W.begin(), true);
  W.R[jpart] = W.R[ipart];  // move j on top of i
  for (int ir=1;ir<nstep;ir++)
  { // slowly move away
    r = ir*dr;
    W.R[jpart][0] = W.R[ipart][0] + r;
    W.update();
    logpsi = Psi.evaluateLog(W);
    eloc   = H.evaluate(W);
    eout << r << " " << eloc << " " << std::real(logpsi) << "\n";
  }
  //MCWalkerConfiguration::iterator Wit(W.begin());
  //for (; Wit != W.end(); Wit++)
  //{
  //  W.loadWalker(**Wit, true);
  //  break;
  //}
}
bool WaveFunctionPlotter::put(xmlNodePtr q)
{
  //set defaults
  nstep = 100;
  ipart = 0;
  jpart = -1;
  // read inputs
  m_param.add(nstep, "nstep");
  m_param.add(ipart, "i");
  m_param.add(jpart, "j");
  m_param.put(q);
  // nested inputs
  xmlNodePtr tcur = q->children;
  while (tcur != NULL)
  {
    std::string cname((const char*)(tcur->name));
    //if (cname)
    //{
    //}
    tcur = tcur->next;
  }
  bool success = putQMCInfo(q);
  return success;
}
} // namespace qmcplusplus
