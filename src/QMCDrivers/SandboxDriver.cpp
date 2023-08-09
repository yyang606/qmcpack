#include "QMCDrivers/SandboxDriver.h"

namespace qmcplusplus
{

SandboxDriver::SandboxDriver(const ProjectData& project_data,
                                   MCWalkerConfiguration& w,
                                   TrialWaveFunction& psi,
                                   QMCHamiltonian& h,
                                   ParticleSetPool& ptclPool,
                                   Communicate* comm)
  : QMCDriver(project_data, w, psi, h, comm, "SandboxDriver"),
  PtclPool(ptclPool),
  ndim(w.getLattice().ndim)
{
}

bool SandboxDriver::run()
{
  plot_obs("moire");
  return true;
}

void SandboxDriver::plot_obs(std::string name)
{
  std::vector<size_t> mesh = {52, 88, 1};
  int nat = W.getTotalNum();
  auto axes = W.getLattice().R;
  //pick the first walker
  Walker_t& awalker = **W.begin();
  W.loadWalker(awalker, false);
  //find observable
  int iobs=-1;
  for (int i = 0; i < H.sizeOfObservables(); i++)
    if (H.getObservableName(i) == name)
      iobs=i;
  if (iobs < 0) throw std::runtime_error(name);
  //march through simulation cell
  std::ofstream fout((name+".dat").c_str());
  for (size_t ix=0;ix<mesh[0];ix++)
  {
    for (size_t iy=0;iy<mesh[1];iy++)
    {
      for (size_t iz=0;iz<mesh[2];iz++)
      {
        PosType r = (RealType)ix*axes.getRow(0)/(RealType)mesh[0];
        r += (RealType)iy*axes.getRow(1)/(RealType)mesh[1];
        r += (RealType)iz*axes.getRow(2)/(RealType)mesh[2];
        W.R = r;
        W.update();
        H.evaluate(W);
        H.auxHevaluate(W);
        auto val = H.getObservable(iobs);
        fout << r << " " << val/nat << std::endl;
      }
    } // y
  } // x
  fout.close();
}

}
