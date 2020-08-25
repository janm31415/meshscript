#include "mm.h"
#include <jtk/concurrency.h>
#include <jtk/file_utils.h>
#include <jtk/qbvh.h>
#include <jtk/geometry.h>

#include <hdf5.h>

#include "hdf5aux.h"
#include "mesh.h"

namespace
  {
  bool read_morphable_model_basel(hid_t group_id, jtk::morphable_model& mm, const char* filename)
    {
    using namespace jtk;   
    herr_t status;

    auto represented_id = H5Gopen(group_id, "representer", H5P_DEFAULT);

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
      mm.triangles[t][2] = data[t + 2 * maxdims[1]];
      }

    status = H5Dclose(cells_id);
    status = H5Sclose(cells_space_id);
    status = H5Gclose(represented_id);

    auto model_id = H5Gopen(group_id, "model", H5P_DEFAULT);

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
      s = std::sqrt(s)*std::sqrt(mm.S.rows() - 1);

    status = H5Dclose(variance_id);
    status = H5Sclose(variance_space_id);

    mm.V.resize(0, 0);

    status = H5Gclose(model_id);
    return true;
    }

  bool read_morphable_model_basel(jtk::morphable_model& shape, jtk::morphable_model& color, const char* filename)
    {
    using namespace jtk;

    hid_t file_id;
    herr_t status;
    file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    auto shape_id = H5Gopen(file_id, "shape", H5P_DEFAULT);
    bool result = read_morphable_model_basel(shape_id, shape, filename);   
    status = H5Gclose(shape_id);

    auto color_id = H5Gopen(file_id, "color", H5P_DEFAULT);
    result &= read_morphable_model_basel(color_id, color, filename);
    status = H5Gclose(color_id);
   
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
    if (!jtk::read_morphable_model_binary(morph.shape, filename.c_str()))
      return false;
    }
  else if (ext == "h5")
    {
    if (!read_morphable_model_basel(morph.shape, morph.color, filename.c_str()))
      return false;
    }
  else
    return false;
  morph.coefficients.resize(morph.shape.S.rows(), 0.f);
  morph.vertices = jtk::get_vertices(morph.shape, morph.coefficients);  
  morph.color_coefficients.resize(morph.color.S.rows(), 0.f);
  morph.vertex_colors = jtk::get_vertices(morph.color, morph.color_coefficients);
  clamp_vertex_colors(morph.vertex_colors);
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
  for (const auto& tria : m.shape.triangles)
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
    return jtk::write_morphable_model_binary(morph.shape, filename.c_str());
    }  
  return false;
  }

void clamp_vertex_colors(std::vector<jtk::vec3<float>>& vertex_colors)
  {
  for (auto& v : vertex_colors)
    {
    for (int j = 0; j < 3; ++j)
      v[j] = v[j] < 0.f ? 0.f : (v[j] > 1.f ? 1.f : v[j]);
    }
  }

void fit_to_mesh(mm& morph, const mesh& m)
  {
  using namespace jtk;
  std::vector<vec3<float>> shape = morph.vertices;
  uint32_t sz = (uint32_t)morph.vertices.size();

  qbvh bvh(m.triangles, m.vertices.data());

  float4x4 cs = m.cs;
  float4x4 cs_inv = invert_orthonormal(m.cs);
  int max_iter = 5;
  for (int iter = 0; iter < max_iter; ++iter)
    {
    std::vector<jtk::vec3<float>> triangle_normals;
    compute_triangle_normals(triangle_normals, morph.vertices.data(), morph.shape.triangles.data(), (uint32_t)morph.shape.triangles.size());

    parallel_for(uint32_t(0), sz, [&](uint32_t i)
      {
      float4 pos(shape[i][0], shape[i][1], shape[i][2], 1.f);
      pos = matrix_vector_multiply(cs_inv, pos);
      float4 dir(triangle_normals[i][0], triangle_normals[i][1], triangle_normals[i][2], 0.f);
      dir = matrix_vector_multiply(cs_inv, dir);

      ray r;
      r.orig = pos;
      r.dir = dir;
      r.t_near = -200.f;
      r.t_far = std::numeric_limits<float>::max();

      uint32_t triangle_id;
      //vec3<float> pt(pos[0], pos[1], pos[2]);
      //auto h = bvh.find_closest_triangle(triangle_id, pt, m.triangles.data(), m.vertices.data());

      auto h = bvh.find_closest_triangle(triangle_id, r, m.triangles.data(), m.vertices.data());

      if (h.found)
        {
        const auto v0 = m.vertices[m.triangles[triangle_id][0]];
        const auto v1 = m.vertices[m.triangles[triangle_id][1]];
        const auto v2 = m.vertices[m.triangles[triangle_id][2]];

        const auto shape_vert = v0 * (1.f - h.u - h.v) + h.u*v1 + h.v*v2;
        float4 sv(shape_vert[0], shape_vert[1], shape_vert[2], 1.f);
        sv = matrix_vector_multiply(cs, sv);
        shape[i][0] = sv[0];
        shape[i][1] = sv[1];
        shape[i][2] = sv[2];
        }
      });
    bool sigma_constraint = iter > max_iter / 2;
    sigma_constraint = false;
    std::vector<float> coeff = fit_shape(shape, morph.shape, sigma_constraint);
    morph.coefficients = coeff;


    morph.vertices = jtk::get_vertices(morph.shape, morph.coefficients);
    morph.vertex_colors = jtk::get_vertices(morph.color, morph.color_coefficients);

    shape = morph.vertices;
    }
  clamp_vertex_colors(morph.vertex_colors);
  /*
  std::vector<float> coeff = fit_shape(shape, morph.shape, false);
  morph.coefficients = coeff;
  
  
  morph.vertices = jtk::get_vertices(morph.shape, morph.coefficients);
  morph.vertex_colors = jtk::get_vertices(morph.color, morph.color_coefficients);
  */
  }

void fit_to_partial_positions(mm& morph, const std::vector<uint32_t>& vertex_indices, const std::vector<jtk::vec3<float>>& vertex_positions)
  {
  if (vertex_indices.size() != vertex_positions.size())
    return;

  using namespace jtk;

  int max_iter = 5;
  for (int iter = 0; iter < max_iter; ++iter)
    {
    std::vector<vec3<float>> shape = morph.vertices;
    for (uint32_t i = 0; i < (uint32_t)vertex_indices.size(); ++i)
      {
      uint32_t idx = vertex_indices[i];
      if (std::isnan(vertex_positions[i][0]))
        continue;
      if (std::isnan(vertex_positions[i][1]))
        continue;
      if (std::isnan(vertex_positions[i][2]))
        continue;
      if (idx < shape.size())
        {
        //shape[idx] = shape[idx]+(vertex_positions[i] - shape[idx])*pow((float)(max_iter-iter), 2.f);
        shape[idx] = vertex_positions[i];
        }
      }
    bool sigma_constraint = iter > max_iter / 2;
    sigma_constraint = false;
    std::vector<float> coeff = fit_shape(shape, morph.shape, sigma_constraint);
    morph.coefficients = coeff;

    morph.vertices = jtk::get_vertices(morph.shape, morph.coefficients);
    //morph.vertex_colors = jtk::get_vertices(morph.color, morph.color_coefficients);
    }
  }

void fit(mm& morph, const mesh& m, const std::vector<uint32_t>& vertex_indices, const std::vector<jtk::vec3<float>>& vertex_positions)
  {
  using namespace jtk;
  std::vector<vec3<float>> shape = morph.vertices;
  std::vector<vec3<float>> color = morph.vertex_colors;
  uint32_t sz = (uint32_t)morph.vertices.size();

  qbvh bvh(m.triangles, m.vertices.data());

  float4x4 cs = m.cs;
  float4x4 cs_inv = invert(m.cs);
  int max_iter = 5;
  for (int iter = 0; iter < max_iter; ++iter)
    {
    std::vector<jtk::vec3<float>> triangle_normals;
    compute_triangle_normals(triangle_normals, morph.vertices.data(), morph.shape.triangles.data(), (uint32_t)morph.shape.triangles.size());

    parallel_for(uint32_t(0), sz, [&](uint32_t i)
      {
      float4 pos(shape[i][0], shape[i][1], shape[i][2], 1.f);
      pos = matrix_vector_multiply(cs_inv, pos);
      float4 dir(triangle_normals[i][0], triangle_normals[i][1], triangle_normals[i][2], 0.f);
      dir = matrix_vector_multiply(cs_inv, dir);

      ray r;
      r.orig = pos;
      r.dir = dir;
      r.t_near = -200.f;
      r.t_far = std::numeric_limits<float>::max();

      uint32_t triangle_id;
      //vec3<float> pt(pos[0], pos[1], pos[2]);
      //auto h = bvh.find_closest_triangle(triangle_id, pt, m.triangles.data(), m.vertices.data());

      auto h = bvh.find_closest_triangle(triangle_id, r, m.triangles.data(), m.vertices.data());

      if (h.found)
        {
        const auto v0 = m.vertices[m.triangles[triangle_id][0]];
        const auto v1 = m.vertices[m.triangles[triangle_id][1]];
        const auto v2 = m.vertices[m.triangles[triangle_id][2]];

        const auto shape_vert = v0 * (1.f - h.u - h.v) + h.u*v1 + h.v*v2;
        float4 sv(shape_vert[0], shape_vert[1], shape_vert[2], 1.f);
        sv = matrix_vector_multiply(cs, sv);
        shape[i][0] = sv[0];
        shape[i][1] = sv[1];
        shape[i][2] = sv[2];

        if (!color.empty() && !m.vertex_colors.empty())
          {
          const auto c0 = m.vertex_colors[m.triangles[triangle_id][0]];
          const auto c1 = m.vertex_colors[m.triangles[triangle_id][1]];
          const auto c2 = m.vertex_colors[m.triangles[triangle_id][2]];
          const auto new_color = c0 * (1.f - h.u - h.v) + h.u*c1 + h.v*c2;
          color[i] = new_color;
          }
        }
      });

    for (uint32_t i = 0; i < (uint32_t)vertex_indices.size(); ++i)
      {
      uint32_t idx = vertex_indices[i];
      if (std::isnan(vertex_positions[i][0]))
        continue;
      if (std::isnan(vertex_positions[i][1]))
        continue;
      if (std::isnan(vertex_positions[i][2]))
        continue;
      if (idx < shape.size())
        shape[idx] = shape[idx] + (vertex_positions[i] - shape[idx])*(float)(max_iter-iter)*(float)(max_iter - iter);
      }

    bool sigma_constraint = iter > max_iter / 2;
    sigma_constraint = false;
    std::vector<float> coeff = fit_shape(shape, morph.shape, sigma_constraint);
    morph.coefficients = coeff;
    morph.color_coefficients = fit_shape(color, morph.color, sigma_constraint);

    morph.vertices = jtk::get_vertices(morph.shape, morph.coefficients);
    morph.vertex_colors = jtk::get_vertices(morph.color, morph.color_coefficients);

    shape = morph.vertices;
    color = morph.vertex_colors;
    }
  clamp_vertex_colors(morph.vertex_colors);
  }