#include "pc.h"
#include "io.h"

#include <jtk/file_utils.h>
#include <jtk/geometry.h>
#include <iostream>

#include <algorithm>

using namespace jtk;

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
    if (p.normals.empty() && p.vertex_colors.empty())
      return jtk::write_ply(filename.c_str(), p.vertices);
    else if (p.vertex_colors.empty())
      return jtk::write_ply(filename.c_str(), p.vertices, p.normals);
    else
      return jtk::write_ply(filename.c_str(), p.vertices, p.normals, p.vertex_colors);
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