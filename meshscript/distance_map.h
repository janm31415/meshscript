#pragma once

#include <jtk/qbvh.h>
#include <jtk/vec.h>
#include <vector>

std::vector<float> distance_map(
  const std::vector<jtk::vec3<float>>& source_vertices,
  const jtk::float4x4& source_cs,
  const std::vector<jtk::vec3<float>>& target_vertices,
  const std::vector<jtk::vec3<uint32_t>>& target_triangles,
  const jtk::float4x4& target_cs,
  bool signed_distance);