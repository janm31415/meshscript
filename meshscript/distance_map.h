#pragma once

#include <jtk/vec.h>
#include <vector>

std::vector<float> distance_map(
  const std::vector<jtk::vec3<float>>& source_vertices,
  const std::vector<jtk::vec3<float>>& target_vertices,
  const std::vector<jtk::vec3<uint32_t>>& target_triangles,
  bool signed_distance);