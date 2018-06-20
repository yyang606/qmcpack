#include <fstream>
#include <sstream>
#include "Message/Communicate.h"
#include "QMCApp/QMCAppBase.h"
#include "Lattice/Uniform3DGridLayout.h"
#include "ParticleIO/ParticleLayoutIO.h"
//#include "LongRange/KContainer.h"
#include "Numerics/OneDimGridBase.h"
#include "einspline/bspline_create.h"
#include "einspline/bspline_eval_d.h"

using namespace qmcplusplus;
using namespace std;

typedef QMCTraits::RealType   RealType;
typedef QMCTraits::IndexType IndexType;
typedef QMCTraits::PosType     PosType;
typedef LinearGrid<RealType>  GridType;


xmlNodePtr find(const char* expression, xmlXPathContextPtr context)
{ // find the first node matching xpath expression in the given context
  OhmmsXPathObject xpath(expression, context);
  if (xpath.size() != 1)
  {
    APP_ABORT("expected 1 " << expression << " found " << xpath.size());
  }
  return xpath[0];
}


vector<vector<RealType>> loadtxt(const string fname)
{ // stack overflow "How to read in a data file of unknown dimensions in C/C++"

  vector<vector<RealType>> mat;
  ifstream fdat;
  fdat.open(fname.c_str(), ios_base::in);

  double buf;
  string line;

  while (!fdat.eof())
  {
    getline(fdat, line);

    // skip comment or empty line
    if (line[0] == '#' || line.empty()) continue;

    // read a row of numbers
    vector<RealType> row;
    stringstream ss(line, ios_base::in);
    while (ss>>buf) row.push_back(buf);

    // add a row to matrix
    mat.push_back(row);
  }

  fdat.close();
  return mat;
}


BCtype_d natural_boundary()
{
  BCtype_d bc;
  bc.lCode = NATURAL;
  bc.rCode = NATURAL;
  bc.lVal = 1.0;
  bc.rVal = 1.0;
  return bc;
}


UBspline_3d_d* create_spline3d(
  GridType gridx,
  GridType gridy,
  GridType gridz,
  vector<RealType> sk)
{ // spline 3D volumetric data sk on regular grid
  Ugrid esgridx, esgridy, esgridz;
  esgridx = gridx.einspline_grid();
  esgridy = gridy.einspline_grid();
  esgridz = gridz.einspline_grid();
  BCtype_d bcx, bcy, bcz;
  bcx = natural_boundary();
  bcy = natural_boundary();
  bcz = natural_boundary();
  UBspline_3d_d* spline3d = create_UBspline_3d_d(
    esgridx, esgridy, esgridz
  , bcx, bcy, bcz, sk.data()
  );
  return spline3d;
}


vector<int> get_index3d(
  const vector<double> gvec,
  const vector<double> gmin,
  const vector<double> dg
)
{
  int ix, iy, iz;
  ix = nearbyint((gvec[0]-gmin[0]-dg[0]/2.)/dg[0]);
  iy = nearbyint((gvec[1]-gmin[1]-dg[1]/2.)/dg[1]);
  iz = nearbyint((gvec[2]-gmin[2]-dg[2]/2.)/dg[2]);
  vector<int> idx3d = {ix, iy, iz};
  return idx3d;
}


int index3d_to_index1d(
  const vector<int> idx3d,
  const vector<int> ng
)
{
  int cidx1d = idx3d[0]*ng[1]*ng[2] + idx3d[1]*ng[2] + idx3d[2];
  return cidx1d;
}


vector<int> index1d_to_index3d(const int idx1d, const vector<int> ng)
{
  vector<int> idx3d(3);
  idx3d[2] = idx1d % ng[2];
  idx3d[1] = (idx1d/ng[2]) % ng[1];
  idx3d[0] = ((idx1d/ng[2])/ng[1]) % ng[0];
  return idx3d;
}


int main(int argc, char **argv)
{
  OHMMS::Controller->initialize(argc, argv);
  if (argc<2) APP_ABORT("usage: " << argv[0] << " axes.xml");
 
  // step 1: read <simulationcell> from input
  Libxml2Document fxml;
  fxml.parse(argv[1]);
  xmlXPathContextPtr doc = fxml.getXPathContext();
  xmlNodePtr sc_node = find("//simulationcell", doc);

  // step 2: construct simulation box
  Uniform3DGridLayout box;
  LatticeParser parser(box);
  parser.put(sc_node);
  box.SetLRCutoffs();  // box.LR_rc, box.LR_kc
  RealType volume = box.Volume;
  app_log() << "  lattice volume = " << volume << endl;

  //// step 3: set up simulation kvectors
  RealType kc = box.LR_kc;
  //KContainer kvecs;
  //kvecs.UpdateKLists(box, kc);
  
  // step 3: set up kvectors
  vector<double> gmin, gmax;
  vector<int> ng;
  xmlNodePtr node;
  node = find("//gmin", doc);
  putContent(gmin, node);
  node = find("//gmax", doc);
  putContent(gmax, node);
  node = find("//ng", doc);
  putContent(ng, node);
  vector<double> dg(3, 0.0);

  app_log() << "using user-defined grid in reciprocal lattice units" << endl;
  app_log() << "input grid: gmin gmax ng dg" << endl;
  for (int idim=0; idim<3; idim++)
  {
    dg[idim] = (gmax[idim]-gmin[idim])/(ng[idim]-1);
    app_log() << gmin[idim] << " " << gmax[idim] << " " << ng[idim] << " " << dg[idim] << endl;
  }
  // initialize kgrid and values
  GridType gridx, gridy, gridz;
  gridx.set(gmin[0], gmax[0], ng[0]);
  gridy.set(gmin[1], gmax[1], ng[1]);
  gridz.set(gmin[2], gmax[2], ng[2]);

  vector<RealType> sk(ng[0]*ng[1]*ng[2], -1);

  // step 4: load data on simulation kgrid
  vector<vector<RealType>> mat = loadtxt("sofk.dat");

  // step 5: transfer data to regular grid
  int nk = mat.size();
  RealType maxval=mat[0][3];
  for (int ik=0; ik<nk; ik++)
  {
    PosType kvec(mat[ik][0], mat[ik][1], mat[ik][2]);
    RealType val = mat[ik][3];
    PosType gvec = box.k_unit(kvec);

    vector<double> mygvec = {gvec[0], gvec[1], gvec[2]};
    vector<int> idx3d = get_index3d(mygvec, gmin, dg);
    int cidx1d = index3d_to_index1d(idx3d, ng);
    sk[cidx1d] = val;
    if (val>maxval) maxval=val;
  }

  // !!!! set unavailable spots
  ofstream ofk;
  ofk.open("missing_kvecs.dat");
  for (int isk=0; isk<sk.size(); isk++)
  {
    if (sk[isk]==-1){
      sk[isk] = 0;

      vector<int> idx3d = index1d_to_index3d(isk, ng);
      PosType gvec(
        gridx[idx3d[0]],
        gridy[idx3d[1]],
        gridz[idx3d[2]]
      );
      PosType kvec = box.k_cart(gvec);
      ofk << kvec << " " << isk << endl;
    }
  }
  ofk.close();
  // !!!! set k=0 to 0
  vector<double> gvec = {0, 0, 0};
  vector<int> idx3d = get_index3d(gvec, gmin, dg);
  int cidx1d = index3d_to_index1d(idx3d, ng);
  sk[cidx1d] = maxval;


  // output S(k) on regular grid for debugging
  ofstream ofs;
  ofs.open("skgrid.dat", ofstream::out);
  for (int ix=0; ix<ng[0]; ix++)
  {
    for (int iy=0; iy<ng[1]; iy++)
    {
      for (int iz=0; iz<ng[2]; iz++)
      {
        int cidx1d = ix*ng[1]*ng[2] + iy*ng[2] + iz;
        RealType val = sk[cidx1d];
        PosType kvec(gridx[ix], gridy[iy], gridz[iz]);
        ofs << kvec << " " << val << endl;
      }
    }
  }
  ofs.close();

  // step 5: create 3D spline
  UBspline_3d_d* spline3d = create_spline3d(gridx, gridy, gridz, sk);

  // dump 3D S(k) spline
  ofs.open("sk3d.dat", ofstream::out);
  int nx = 16;
  RealType kmin, kmax, dk;
  kmin = -0.8;
  kmax =  0.8;
  dk = (kmax-kmin)/(nx-1);

  for (int ix=0; ix<nx; ix++)
  {
    for (int iy=0; iy<nx; iy++)
    {
      for (int iz=0; iz<nx; iz++)
      {
        PosType kvec(
          0.5*dk+kmin+ix*dk,
          0.5*dk+kmin+iy*dk,
          0.5*dk+kmin+iz*dk
        );
        PosType gvec = box.k_unit(kvec);
        RealType val;
        eval_UBspline_3d_d(spline3d, gvec[0], gvec[1], gvec[2], &val);
        ofs << kvec << " " << val << endl;
      }
    }
  }
  ofs.close();
  
  OHMMS::Controller->finalize();
  return 0;
}
