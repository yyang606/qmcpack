//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file observable_helper.h
 *@brief Declaration of observable_helper and other helper class for observables
 */
#ifndef QMCPLUSPLUS_OBSERVABLE_TRAITS_H
#define QMCPLUSPLUS_OBSERVABLE_TRAITS_H

#include "OhmmsData/HDFAttribIO.h"
#include "Numerics/HDFNumericAttrib.h"
#include "Numerics/HDFSTLAttrib.h"
#include "OhmmsData/HDFStringAttrib.h"

namespace qmcplusplus
{
/** define observable_helper
 *
 * This is a helper class to manage a hdf5 dagroup for each collectable.
 * The data handled by an Estimator should be presented by dense N-dim array in C.
 * /observables/title/value
 */
struct observable_helper
{
  typedef QMCTraits::FullPrecRealType value_type;
  ///this can be handled by template argument
  //enum {type_id=H5T_NATIVE_DOUBLE};
  ///starting index
  hsize_t lower_bound = 0;
  ///id of this observable
  hid_t data_id = -1;
  ///dataspace for value
  hid_t space1_id = -1;
  ///id of the value dataset
  hid_t value1_id = -1;
  ///my dimensions
  std::vector<hsize_t> mydims;
  ///maximum dimensions
  std::vector<hsize_t> maxdims;
  ///current dimensions while extending data
  std::vector<hsize_t> curdims;
  ///offsets
  std::vector<hsize_t> offsets;
  ///name of this observable
  std::string group_name;
  ///true: hdf5 handlers are valid via open, false: otherwise
  bool isopened;

  ///default constructor
  observable_helper(const std::string& title = "dummy")
      : data_id(-1), space1_id(-1), value1_id(-1), group_name(title), isopened(false)
  {}

  observable_helper(const observable_helper&) = delete;
  observable_helper& operator=(const observable_helper&) = delete;

  /**
   * Move constructor. Must properly transfer ownership of resources. Doing copies as these are "cheap" elements
   * @param in input object to be moved to this
   */
  observable_helper(observable_helper&& in) noexcept
      : lower_bound(in.lower_bound),
        data_id(in.data_id),
        space1_id(in.space1_id),
        value1_id(in.value1_id),
        mydims(in.mydims),
        maxdims(in.maxdims),
        curdims(in.curdims),
        offsets(in.offsets),
        group_name(in.group_name),
        isopened(in.isopened)
  {
    in.isopened = false;
  }

  observable_helper& operator=(observable_helper&& in)
  {
    if (this != &in)
    {
      *this       = std::move(in);
      in.isopened = false;
    }
    return *this;
  }

  ///close resources
  ~observable_helper()
  {
    if (isopened)
    {
      close();
    }
  }

  /** set the shape of this observable
   * @param dims dimensions
   * @param first starting index
   */
  void set_dimensions(std::vector<int>& dims, int first)
  {
    //rank is increased
    hsize_t rank = dims.size() + 1;
    mydims.resize(rank, 1);
    copy(dims.begin(), dims.end(), mydims.begin() + 1);
    maxdims = mydims;
    curdims = mydims;
    offsets.resize(rank, 0);
    maxdims[0]  = H5S_UNLIMITED;
    lower_bound = first;
  }

  /** open a h5 group of this observable
   *
   * Create a group for an observable and dataspace
   */
  inline void open(hid_t grp_id)
  {
    data_id      = H5Gcreate(grp_id, group_name.c_str(), 0);
    hsize_t rank = mydims.size();
    if (rank)
    {
      //create empty data to write something
      hsize_t nd = 1;
      for (int i = 1; i < rank; ++i)
        nd *= mydims[i];
      std::vector<value_type> zeros(nd, 0.0);
      hid_t p = H5Pcreate(H5P_DATASET_CREATE);
      H5Pset_chunk(p, rank, &mydims[0]);
      space1_id      = H5Screate_simple(rank, &mydims[0], &maxdims[0]);
      value1_id      = H5Dcreate(data_id, "value", H5T_NATIVE_DOUBLE, space1_id, p);
      hid_t memspace = H5Screate_simple(rank, &mydims[0], NULL);
      herr_t ret     = H5Dwrite(value1_id, H5T_NATIVE_DOUBLE, memspace, space1_id, H5P_DEFAULT, &zeros[0]);
      H5Sclose(memspace);
      H5Pclose(p);
    }
    isopened = true;
  }

  void close()
  {
    if (isopened)
    {
      if (space1_id > -1)
      {
        H5Sclose(space1_id);
        space1_id = -1;
      }
      if (data_id > -1)
      {
        H5Gclose(data_id);
        data_id = -1;
      }
      isopened = false;
    }
  }

  /** add named property to describe the collectable this helper class handles
   * @param p any intrinsic datatype including vector, basic containers
   * @param pname property
   */
  template<typename T>
  inline void addProperty(T& p, const std::string& pname)
  {
    HDFAttribIO<T> a(p);
    a.write(data_id, pname.c_str());
  }

  inline void addProperty(float& p, const std::string& pname)
  {
    double p_DP(p);
    HDFAttribIO<double> a(p_DP);
    a.write(data_id, pname.c_str());
  }

  inline void addProperty(Tensor<float, OHMMS_DIM>& p, const std::string& pname)
  {
    Tensor<double, OHMMS_DIM> p_DP;
    p_DP = p;
    HDFAttribIO<Tensor<double, OHMMS_DIM>> a(p_DP);
    a.write(data_id, pname.c_str());
  }

  inline void addProperty(Matrix<float>& p, const std::string& pname)
  {
    Matrix<double> p_DP;
    p_DP = p;
    HDFAttribIO<Matrix<double>> a(p_DP);
    a.write(data_id, pname.c_str());
  }

  inline void addProperty(TinyVector<float, OHMMS_DIM>& p, const std::string& pname)
  {
    TinyVector<double, OHMMS_DIM> p_DP(p);
    HDFAttribIO<TinyVector<double, OHMMS_DIM>> a(p_DP);
    a.write(data_id, pname.c_str());
  }

  inline void addProperty(std::vector<float>& p, const std::string& pname)
  {
    std::vector<double> p_DP;
    p_DP.assign(p.begin(), p.end());
    HDFAttribIO<std::vector<double>> a(p_DP);
    a.write(data_id, pname.c_str());
  }

  inline void addProperty(std::vector<TinyVector<float, OHMMS_DIM>>& p, const std::string& pname)
  {
    std::vector<TinyVector<double, OHMMS_DIM>> p_DP;
    p_DP.assign(p.begin(), p.end());
    HDFAttribIO<std::vector<TinyVector<double, OHMMS_DIM>>> a(p_DP);
    a.write(data_id, pname.c_str());
  }

  inline void write(const value_type* first_v, const value_type* first_vv)
  {
    hsize_t rank = mydims.size();
    if (rank)
    {
      H5Sset_extent_simple(space1_id, rank, &curdims[0], &maxdims[0]);
      H5Sselect_hyperslab(space1_id, H5S_SELECT_SET, &offsets[0], NULL, &mydims[0], NULL);
      H5Dextend(value1_id, &curdims[0]);
      hid_t memspace = H5Screate_simple(rank, &mydims[0], NULL);
      herr_t ret     = H5Dwrite(value1_id, H5T_NATIVE_DOUBLE, memspace, space1_id, H5P_DEFAULT, first_v + lower_bound);
      H5Sclose(memspace);
      curdims[0]++;
      offsets[0]++;
    }
  }
};
} // namespace qmcplusplus
#endif
