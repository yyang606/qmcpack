#include "WaveFunctionPlotter.h"

namespace qmcplusplus
{

using std::endl;

double read_column(ParticleSet& P, std::string name)
{
  double value = -9999999;
  bool found = false;
  for (int j=0;j<P.PropertyList.size();j++)
  {
    if (name == P.PropertyList.Names[j])
    {
      value = P.PropertyList.Values[j];
      found = true;
      break;
    }
  }
  if (not found) throw value;
  return value;
}

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
  app_log() << npart << " particles" << endl;
  // load first walker
  W.loadWalker(**W.begin(), true);
  // put particle at 0
  W.R[0] = 0;
  for (int i=0;i<npart;i++)
  {
    app_log() << W.R[i] << endl;
  }
  //H.evaluate(W);
  H.auxHevaluate(W);
  app_log() << read_column(W, "moire") << endl;

  double dx=0.25, dy=0.25, x, y, v;
  int nx=32, ny=32;
  //int nx=2, ny=2;
  std::ofstream fout("vpot.dat");
  for (int ix=0;ix<nx;ix++)
  {
    for (int iy=0;iy<ny;iy++)
    {
      x = ix*dx;
      y = iy*dy;
      app_log() << x << " " << y << endl;
      W.R[0] = {x, y, 0.0};
      W.update();
      H.auxHevaluate(W);
      v = read_column(W, "moire");
      fout << x << " " << y << " " << v << endl;
    }
  }
  fout.close();

  //for (int j=0;j<W.PropertyList.size();j++)
  //{
  //  app_log() << W.PropertyList.Names[j] << " " << W.PropertyList.Values[j] << endl;
  //}
  /*
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
  */
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
