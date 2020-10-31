#include "pc.h"
#include "io.h"

#include <jtk/file_utils.h>
#include <jtk/fitting.h>
#include <jtk/geometry.h>
#include <jtk/mat.h>
#include <iostream>

#include <algorithm>

#include <jtk/point_tree.h>

#include <jtk/containers.h>

#include <trico/trico/trico.h>

using namespace jtk;

namespace
  {

  bool write_trc(const pc& p, const std::string& filename)
    {
    void* arch = trico_open_archive_for_writing(1024 * 1024);
    if (!trico_write_vertices(arch, (float*)p.vertices.data(), (uint32_t)p.vertices.size()))
      {
      std::cout << "Something went wrong when writing the vertices\n";
      return false;
      }
    if (!p.vertex_colors.empty() && !trico_write_vertex_colors(arch, (uint32_t*)p.vertex_colors.data(), (uint32_t)p.vertex_colors.size()))
      {
      std::cout << "Something went wrong when writing the vertex colors\n";
      return false;
      }
    if (!p.normals.empty() && !trico_write_vertex_normals(arch, (float*)p.normals.data(), (uint32_t)p.normals.size()))
      {
      std::cout << "Something went wrong when writing the normals\n";
      return false;
      }

    FILE* f = fopen(filename.c_str(), "wb");
    if (!f)
      {
      std::cout << "Cannot write to file " << filename << std::endl;
      return false;
      }

    fwrite((const void*)trico_get_buffer_pointer(arch), trico_get_size(arch), 1, f);
    fclose(f);

    trico_close_archive(arch);

    return true;
    }

  bool read_trc(pc& p, const std::string& filename)
    {
    long long size = jtk::file_size(filename.c_str());
    if (size < 0)
      {
      std::cout << "There was an error reading file " << filename << std::endl;
      return false;
      }
    FILE* f = fopen(filename.c_str(), "rb");
    if (!f)
      {
      std::cout << "Cannot open file: " << filename << std::endl;
      return false;
      }
    char* buffer = (char*)malloc(size);
    long long fl = (long long)fread(buffer, 1, size, f);
    if (fl != size)
      {
      std::cout << "There was an error reading file " << filename << std::endl;
      fclose(f);
      return false;
      }
    fclose(f);

    void* arch = trico_open_archive_for_reading((const uint8_t*)buffer, size);
    if (!arch)
      {
      std::cout << "The input file " << filename << " is not a trico archive." << std::endl;
      return false;
      }

    enum trico_stream_type st = trico_get_next_stream_type(arch);
    while (!st == trico_empty)
      {
      switch (st)
        {
        case trico_vertex_float_stream:
        {
        p.vertices.resize(trico_get_number_of_vertices(arch));
        float* vertices = (float*)p.vertices.data();
        if (!trico_read_vertices(arch, &vertices))
          {
          std::cout << "Something went wrong reading the vertices" << std::endl;
          trico_close_archive(arch);
          free(buffer);
          return false;
          }
        break;
        }        
        case trico_vertex_color_stream:
        {       
        p.vertex_colors.resize(trico_get_number_of_colors(arch));
        uint32_t* vertex_colors = (uint32_t*)p.vertex_colors.data();
        if (!trico_read_vertex_colors(arch, &vertex_colors))
          {
          std::cout << "Something went wrong reading the vertex colors" << std::endl;
          trico_close_archive(arch);
          free(buffer);
          return false;
          }
        break;
        }
        case trico_vertex_normal_float_stream:
        {
        p.normals.resize(trico_get_number_of_normals(arch));
        float* normals = (float*)p.normals.data();
        if (!trico_read_vertex_normals(arch, &normals))
          {
          std::cout << "Something went wrong reading the normals" << std::endl;
          trico_close_archive(arch);
          free(buffer);
          return false;
          }        
        break;
        }
        default:
        {
        trico_skip_next_stream(arch);
        break;
        }
        }
      st = trico_get_next_stream_type(arch);
      }

    trico_close_archive(arch);
    free(buffer);

    return true;
    }
  }

bool read_from_file(pc& point_cloud, const std::string& filename)
  {
  std::string ext = jtk::get_extension(filename);
  if (ext.empty())
    return false;
  std::transform(ext.begin(), ext.end(), ext.begin(), [](char ch) {return (char)::tolower(ch); });
  if (ext == "ply")
    {
    std::vector<jtk::vec3<uint32_t>> triangles;
    std::vector<jtk::vec3<jtk::vec2<float>>> uv;
    if (!read_ply(filename.c_str(), point_cloud.vertices, point_cloud.normals, point_cloud.vertex_colors, triangles, uv))
      return false;
    if (point_cloud.vertices.empty())
      return false;
    }
  else if (ext == "trc")
    {
    if (!read_trc(point_cloud, filename))
      return false;
    }
  point_cloud.cs = get_identity();
  point_cloud.visible = true;
  return true;
  }

bool vertices_to_csv(const pc& m, const std::string& filename)
  {
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

bool write_to_file(const pc& p, const std::string& filename)
  {
  std::string ext = jtk::get_extension(filename);
  if (ext.empty())
    return false;
  if (ext == "ply")
    {
    /*
    if (p.normals.empty() && p.vertex_colors.empty())
      return jtk::write_ply(filename.c_str(), p.vertices);
    else if (p.vertex_colors.empty())
      return jtk::write_ply(filename.c_str(), p.vertices, p.normals);
    else
      return jtk::write_ply(filename.c_str(), p.vertices, p.normals, p.vertex_colors);
    */
    std::vector<jtk::vec3<uint32_t>> triangles;
    std::vector<jtk::vec3<jtk::vec2<float>>> uv;
    return write_ply(filename.c_str(), p.vertices, p.normals, p.vertex_colors, triangles, uv);
    }
  else if (ext == "trc")
    {
    return write_trc(p, filename);
    }
  return false;
  }

void info(const pc& p)
  {
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "POINTCLOUD" << std::endl;
  std::cout << "Vertices: " << p.vertices.size() << std::endl;
  std::cout << "Coordinate system: " << std::endl;
  for (int i = 0; i < 4; ++i)
    {
    for (int j = 0; j < 4; ++j)
      {
      std::cout << p.cs[i + 4 * j] << " ";
      }
    std::cout << std::endl;
    }
  std::cout << "Vertex normals: " << (p.normals.empty() ? "No" : "Yes") << std::endl;
  std::cout << "Vertex colors: " << (p.vertex_colors.empty() ? "No" : "Yes") << std::endl;
  std::cout << "Visible: " << (p.visible ? "Yes" : "No") << std::endl;
  std::cout << "---------------------------------------" << std::endl;
  }

void cs_apply(pc& p)
  {
  for (auto& v : p.vertices)
    {
    jtk::float4 V(v[0], v[1], v[2], 1.f);
    V = jtk::matrix_vector_multiply(p.cs, V);
    v[0] = V[0];
    v[1] = V[1];
    v[2] = V[2];
    }
  p.cs = jtk::get_identity();
  }

namespace
  {
  struct tree_point
    {
    jtk::vec3<float> pt;
    uint32_t idx;
    float operator [](int i) const
      {
      return pt[i];
      }
    float& operator [](int i)
      {
      return pt[i];
      }
    };

  struct point_tree_traits
    {
    typedef float value_type;
    enum { dimension = 3 };
    typedef tree_point point_type;
    };

  uint64_t make_edge(uint32_t v0, uint32_t v1)
    {
    uint64_t e = v0;
    e <<= 32;
    e |= (uint64_t)v1;
    return e;
    }

  uint32_t edge_to_v0(uint64_t e)
    {
    e >>= 32;
    return (uint32_t)(e & 0xffffffff);
    }

  uint32_t edge_to_v1(uint64_t e)
    {
    return (uint32_t)(e & 0xffffffff);
    }
  }

std::vector<jtk::vec3<float>> estimate_normals(const pc& p, uint32_t k)
  {
  std::vector<jtk::vec3<float>> normals;
  normals.reserve(p.vertices.size());

  point_tree<point_tree_traits> tree;
  std::vector<tree_point> vert;
  vert.reserve(p.vertices.size());
  for (uint32_t v = 0; v < p.vertices.size(); ++v)
    {
    tree_point pt;
    pt.pt = p.vertices[v];
    pt.idx = v;
    vert.push_back(pt);
    }
  tree.efficient_build_tree(vert.begin(), vert.end());

  for (uint32_t v = 0; v < p.vertices.size(); ++v)
    {
    tree_point tp;
    tp.pt = p.vertices[v];;
    std::vector < tree_point > pts = tree.find_k_nearest((int)k, tp);
    jtk::vec3<float> origin, normal;
    float eig;
    std::vector<jtk::vec3<float>> raw_pts;
    raw_pts.reserve(pts.size());
    for (auto p : pts)
      raw_pts.push_back(p.pt);
    jtk::fit_plane(origin, normal, eig, raw_pts);
    normals.push_back(normal);
    }

  std::vector<bool> treated(p.vertices.size(), false);

  jtk::hashed_heap<uint64_t, float> heap;

  uint32_t v = 0;
  while (true)
    {
    while (v < p.vertices.size() && treated[v])
      ++v;
    if (v == p.vertices.size())
      break;

    treated[v] = true;

    tree_point tp;
    tp.pt = p.vertices[v];
    std::vector < tree_point > pts = tree.find_k_nearest((int)k, tp);
    for (const auto& pt : pts)
      {
      if (pt.idx != v && !treated[pt.idx])
        {
        float score = fabs(dot(normals[v], normals[pt.idx]));
        heap.push(std::pair<uint64_t, float>(make_edge(v, pt.idx), score));
        }
      }

    while (!heap.empty())
      {
      auto top_element = heap.top();
      heap.pop();
      uint32_t v0 = edge_to_v0(top_element.first);
      uint32_t v1 = edge_to_v1(top_element.first);
      if (!treated[v1])
        {
        treated[v1] = true;
        if (dot(normals[v0], normals[v1]) < 0.f)
          normals[v1] = -normals[v1];

        tp.pt = p.vertices[v1];
        std::vector < tree_point > pts = tree.find_k_nearest((int)k, tp);
        for (const auto& pt : pts)
          {
          if (pt.idx != v1 && !treated[pt.idx])
            {
            float score = fabs(dot(normals[v1], normals[pt.idx]));
            heap.push(std::pair<uint64_t, float>(make_edge(v1, pt.idx), score));
            }
          }

        }
      }
    }

  return normals;
  }