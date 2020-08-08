#pragma once

#include <jtk/qbvh.h>
#include <jtk/vec.h>

#include <stdint.h>

#include <memory>
#include <string>

struct pc
  {
  std::vector<jtk::vec3<float>> vertices;  
  std::vector<jtk::vec3<float>> normals;
  std::vector<uint32_t> vertex_colors;  
  jtk::float4x4 cs;
  bool visible;
  };

bool read_from_file(pc& point_cloud, const std::string& filename);

bool vertices_to_csv(const pc& m, const std::string& filename);