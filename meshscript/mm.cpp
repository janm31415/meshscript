#include "mm.h"
#include <jtk/file_utils.h>
#include <jtk/qbvh.h>

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