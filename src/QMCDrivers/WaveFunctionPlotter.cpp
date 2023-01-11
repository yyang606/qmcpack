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
  if (not found) throw name;
  return value;
}

WaveFunctionPlotter::WaveFunctionPlotter(const ProjectData& project_data,
                                         MCWalkerConfiguration& w,
                                         TrialWaveFunction& psi,
                                         QMCHamiltonian& h,
                                         Communicate* comm)
    : QMCDriver(project_data, w, psi, h, comm, "WaveFunctionPlotter")
{
}

void WaveFunctionPlotter::plotExternalPotential()
{
  std::ofstream fout(extpot + ".dat");
  int npart = W.getTotalNum();
  auto cell = W.getLattice();
  PosType r;
  std::vector<int> mesh = {135, 135};
  for (int i=0;i<mesh[0];i++)
  {
    for (int j=0;j<mesh[1];j++)
    {
      r = (double)i/mesh[0] * cell.Rv[0] + (double)j/mesh[1] * cell.Rv[1];
      for (int iat=0;iat<npart;iat++) W.R[iat] = r;
      W.update();
      H.evaluate(W);
      H.auxHevaluate(W);
      auto v = read_column(W, extpot);
      fout << r << " " << v/npart << endl;
    }
  }
  fout.close();
}

bool WaveFunctionPlotter::run()
{
  app_log() << "Starting WF plotter\n";
  // load first walker
  W.loadWalker(**W.begin(), true);
  // do something
  if (extpot != "none") plotExternalPotential();
  return true;
}

bool WaveFunctionPlotter::put(xmlNodePtr q)
{
  //set defaults
  extpot = "none";
  // read inputs
  m_param.add(extpot, "extpot");
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
