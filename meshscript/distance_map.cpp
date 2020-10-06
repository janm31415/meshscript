#include "distance_map.h"

#include <jtk/qbvh.h>
#include <jtk/geometry.h>
#include <jtk/concurrency.h>

using namespace jtk;


std::vector<float> distance_map(
  const std::vector<vec3<float>>& source_vertices,
  const std::vector<vec3<float>>& target_vertices,
  const std::vector<vec3<uint32_t>>& target_triangles,
  bool signed_distance)
  {
  std::vector<float> out(source_vertices.size(), std::numeric_limits<float>::infinity());
  std::vector<vec3<float>> triangle_normals;
  compute_triangle_normals(triangle_normals, target_vertices.data(), target_triangles.data(), (uint32_t)target_triangles.size());

  qbvh tree(target_triangles, target_vertices.data());

  uint32_t sz = (uint32_t)source_vertices.size();

  jtk::parallel_for((uint32_t)0, sz, [&](uint32_t i)
    {
    uint32_t triangle_id;

    auto h = tree.find_closest_triangle(triangle_id, source_vertices[i], target_triangles.data(), target_vertices.data());
    if (signed_distance)
      {
      auto n = triangle_normals[triangle_id];
      const vec3<float> V0 = target_vertices[target_triangles[triangle_id][0]];
      const vec3<float> V1 = target_vertices[target_triangles[triangle_id][1]];
      const vec3<float> V2 = target_vertices[target_triangles[triangle_id][2]];
      const auto pos = V0 * (1.f - h.u - h.v) + h.u*V1 + h.v*V2;
      auto v = (source_vertices[i] - pos);
      float sgn = dot(v, n);
      if (sgn < 0.f)
        h.distance = -h.distance;
      }
    out[i] = h.distance >= 0 ? std::sqrt(h.distance) : -std::sqrt(-h.distance);
    });

  return out;
  }