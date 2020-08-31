#include "mesh.h"
#include "io.h"

#include <jtk/geometry.h>

#include <algorithm>

#include <jtk/file_utils.h>

#include <stb_image.h>

using namespace jtk;

void compute_bb(vec3<float>& min, vec3<float>& max, uint32_t nr_of_vertices, const vec3<float>* vertices)
  {
  if (nr_of_vertices == 0)
    return;
  min[0] = vertices[0][0];
  min[1] = vertices[0][1];
  min[2] = vertices[0][2];
  max[0] = min[0];
  max[1] = min[1];
  max[2] = min[2];
  for (uint32_t i = 1; i < nr_of_vertices; ++i)
    {
    min[0] = std::min<float>(min[0], vertices[i][0]);
    min[1] = std::min<float>(min[1], vertices[i][1]);
    min[2] = std::min<float>(min[2], vertices[i][2]);
    max[0] = std::max<float>(max[0], vertices[i][0]);
    max[1] = std::max<float>(max[1], vertices[i][1]);
    max[2] = std::max<float>(max[2], vertices[i][2]);
    }
  }  

bool read_from_file(mesh& m, const std::string& filename)
  {  
  std::string ext = jtk::get_extension(filename);
  if (ext.empty())
    return false;
  std::transform(ext.begin(), ext.end(), ext.begin(), [](char ch) {return (char)::tolower(ch); });
  if (ext == "stl")
    {
    if (!read_stl(m.vertices, m.triangles, filename.c_str()))
      return false;
    }
  else if (ext == "ply")
    {
    std::vector<uint32_t> vertex_colors;
    if (!read_ply(filename.c_str(), m.vertices, m.triangles, vertex_colors))
      return false;
    if (!vertex_colors.empty())
      {
      m.vertex_colors.clear();
      m.vertex_colors.reserve(vertex_colors.size());
      for (uint32_t clr : vertex_colors)
        {
        uint32_t red = clr & 255;
        uint32_t green = (clr >> 8) & 255;
        uint32_t blue = (clr >> 16) & 255;
        m.vertex_colors.emplace_back((float)red / 255.f, (float)green / 255.f, (float)blue / 255.f);
        }
      }
    }
  else if (ext == "off")
    {
    if (!read_off(m.vertices, m.triangles, filename.c_str()))
      return false;
    }
  else if (ext == "obj")
    {
    std::string mtl_filename;
    if (!read_obj(mtl_filename, m.vertices, m.triangles, m.uv_coordinates, filename.c_str()))
      return false;
    if (!mtl_filename.empty())
      {      
      if (!file_exists(mtl_filename))
        {
        mtl_filename = get_folder(filename) + "/" + mtl_filename;
        }
      std::string texture_filename;
      if (read_texture_filename_from_mtl(texture_filename, mtl_filename.c_str()))
        {
        if (!file_exists(texture_filename))
          {
          texture_filename = get_folder(mtl_filename) + "/" + texture_filename;
          }
        int w, h, nr_of_channels;
        unsigned char* im = stbi_load(texture_filename.c_str(), &w, &h, &nr_of_channels, 4);
        if (im)
          {
          m.texture = jtk::span_to_image(w, h, w, (const uint32_t*)im);
          stbi_image_free(im);
          }
        }
      }
    }
  else
    return false;
  m.cs = get_identity();
  m.visible = true;
  return true;
  }

bool vertices_to_csv(const mesh& m, const std::string& filename)
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

bool triangles_to_csv(const mesh& m, const std::string& filename)
  {
  std::vector<std::vector<std::string>> data;
  for (const auto& tria : m.triangles)
    {
    std::vector<std::string> line;
    line.push_back(std::to_string(tria[0]));
    line.push_back(std::to_string(tria[1]));
    line.push_back(std::to_string(tria[2]));
    data.push_back(line);
    }
  return csv_write(data, filename.c_str(), ",");
  }


bool write_to_file(const mesh& m, const std::string& filename)
  {
  std::string ext = jtk::get_extension(filename);
  if (ext.empty())
    return false;
  std::transform(ext.begin(), ext.end(), ext.begin(), [](char ch) {return (char)::tolower(ch); });
  if (ext == "stl")
    {
    return jtk::write_stl(m.vertices.data(), (uint32_t)m.triangles.size(), m.triangles.data(), nullptr, nullptr, filename.c_str());
    }
  else if (ext == "ply")
    {
    if (m.vertex_colors.empty())
      return jtk::write_ply(filename.c_str(), m.vertices, m.triangles);
    else
      {
      std::vector<uint32_t> colors;
      colors.reserve(m.vertex_colors.size());
      for (auto& clr : m.vertex_colors)
        {
        uint32_t r = (uint32_t)(clr[0] * 255.f);
        uint32_t g = (uint32_t)(clr[1] * 255.f);
        uint32_t b = (uint32_t)(clr[2] * 255.f);
        r = r > 255 ? 255 : r;
        g = g > 255 ? 255 : g;
        b = b > 255 ? 255 : b;
        colors.push_back(0xff000000 | (b << 16) | (g << 8) | r);
        }
      return jtk::write_ply(filename.c_str(), m.vertices, colors, m.triangles);
      }
    }
  return false;
  }