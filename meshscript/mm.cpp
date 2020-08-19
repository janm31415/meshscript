#include "mm.h"
#include <jtk/file_utils.h>
#include <jtk/qbvh.h>

#include <hdf5.h>

#include "hdf5aux.h"

namespace
  {
  bool read_morphable_model_basel(jtk::morphable_model& mm, const char* filename)
    {
    using namespace jtk;

    hid_t file_id;
    herr_t status;
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    auto shape_id = H5Gopen(file_id, "shape", H5P_DEFAULT);

    auto represented_id = H5Gopen(shape_id, "representer", H5P_DEFAULT);

    auto cells_id = H5Dopen2(represented_id, "cells", H5P_DEFAULT);
    auto cells_space_id = H5Dget_space(cells_id);

    hsize_t dims[2], maxdims[2];

    status = H5Sget_simple_extent_dims(cells_space_id, dims, maxdims);

    assert(maxdims[0] == 3);
    assert(status == 2);

    std::vector<uint32_t> data(maxdims[0] * maxdims[1]);

    mm.triangles.resize(maxdims[1]);
    status = H5Dread(cells_id, get_hdf5_type<uint32_t>::get_mem_type_id(), cells_space_id, H5S_ALL, H5P_DEFAULT, (void*)data.data());

    for (uint32_t t = 0; t < maxdims[1]; ++t)
      {
      mm.triangles[t][0] = data[t];
      mm.triangles[t][1] = data[t + maxdims[1]];
      mm.triangles[t][2] = data[t + 2*maxdims[1]];
      }

    status = H5Dclose(cells_id);
    status = H5Sclose(cells_space_id);
    status = H5Gclose(represented_id);

    auto model_id = H5Gopen(shape_id, "model", H5P_DEFAULT);

    auto mean_id = H5Dopen(model_id, "mean", H5P_DEFAULT);
    auto mean_space_id = H5Dget_space(mean_id);

    status = H5Sget_simple_extent_dims(mean_space_id, dims, maxdims);
    assert(status == 1);
    mm.average.resize(maxdims[0], 1);
    status = H5Dread(mean_id, get_hdf5_type<float>::get_mem_type_id(), mean_space_id, H5S_ALL, H5P_DEFAULT, (void*)mm.average.data());

    status = H5Dclose(mean_id);
    status = H5Sclose(mean_space_id);

    auto U_id = H5Dopen2(model_id, "pcaBasis", H5P_DEFAULT);
    auto U_space_id = H5Dget_space(U_id);

    status = H5Sget_simple_extent_dims(U_space_id, dims, maxdims);
    assert(status == 2);
    mm.U.resize(maxdims[0], maxdims[1]);
    status = H5Dread(U_id, get_hdf5_type<float>::get_mem_type_id(), U_space_id, H5S_ALL, H5P_DEFAULT, (void*)mm.U.data());

    status = H5Dclose(U_id);
    status = H5Sclose(U_space_id);

    auto variance_id = H5Dopen(model_id, "pcaVariance", H5P_DEFAULT);
    auto variance_space_id = H5Dget_space(variance_id);

    status = H5Sget_simple_extent_dims(variance_space_id, dims, maxdims);
    assert(status == 1);
    mm.S.resize(maxdims[0], 1);
    status = H5Dread(variance_id, get_hdf5_type<float>::get_mem_type_id(), variance_space_id, H5S_ALL, H5P_DEFAULT, (void*)mm.S.data());

    for (auto& s : mm.S)
      s = std::sqrt(s)*std::sqrt(mm.S.rows()-1);

    status = H5Dclose(variance_id);
    status = H5Sclose(variance_space_id);

    mm.V.resize(0, 0);

    status = H5Gclose(model_id);
    status = H5Gclose(shape_id);
    status = H5Fclose(file_id);
    return true;
    }
  }

bool read_from_file(mm& morph, const std::string& filename)
  {
  std::string ext = jtk::get_extension(filename);
  if (ext.empty())
    return false;
  std::transform(ext.begin(), ext.end(), ext.begin(), [](char ch) {return (char)::tolower(ch); });
  if (ext == "ssm")
    {
    if (!jtk::read_morphable_model_binary(morph.m, filename.c_str()))
      return false;
    }
  else if (ext == "h5")
    {
    if (!read_morphable_model_basel(morph.m, filename.c_str()))
      return false;
    }
  else
    return false;
  morph.coefficients.resize(morph.m.S.rows(), 0.f);
  morph.vertices = jtk::get_vertices(morph.m, morph.coefficients);  
  morph.cs = jtk::get_identity();
  morph.visible = true;
  return true;
  }


bool vertices_to_csv(const mm& m, const std::string& filename)
  {
  using namespace jtk;
  std::vector<std::vector<std::string>> data;
  for (const auto& vertex : m.vertices)
    {
    std::vector<std::string> line;
    line.push_back(std::to_string(vertex[0]));
    line.push_back(std::to_string(vertex[1]));
    line.push_back(std::to_string(vertex[2]));
    data.push_back(line);
    }
  return csv_write(data, filename.c_str(), ",");
  }

bool triangles_to_csv(const mm& m, const std::string& filename)
  {
  using namespace jtk;
  std::vector<std::vector<std::string>> data;
  for (const auto& tria : m.m.triangles)
    {
    std::vector<std::string> line;
    line.push_back(std::to_string(tria[0]));
    line.push_back(std::to_string(tria[1]));
    line.push_back(std::to_string(tria[2]));
    data.push_back(line);
    }
  return csv_write(data, filename.c_str(), ",");
  }

bool write_to_file(const mm& morph, const std::string& filename)
  {
  std::string ext = jtk::get_extension(filename);
  if (ext.empty())
    return false;
  std::transform(ext.begin(), ext.end(), ext.begin(), [](char ch) {return (char)::tolower(ch); });
  if (ext == "ssm")
    {
    return jtk::write_morphable_model_binary(morph.m, filename.c_str());
    }  
  return false;
  }