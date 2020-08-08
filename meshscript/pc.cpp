#include "pc.h"
#include "io.h"

#include <jtk/file_utils.h>

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
    if (!read_ply(filename.c_str(), point_cloud.vertices, point_cloud.normals, point_cloud.vertex_colors))
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