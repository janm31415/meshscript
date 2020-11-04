#pragma once

#include <vector>
#include <jtk/vec.h>

struct icosahedron
  {
  icosahedron()
    {
    using namespace jtk;
    vertices.resize(12);
    const float t = (1.f + std::sqrt(5.f)) / 2.f;
    const float denom = std::sqrt(1.f + t * t);
    vertices[0] = vec3<float>(t, 1.0, 0.0) / denom;
    vertices[1] = vec3<float>(-t, 1.0, 0.0) / denom;
    vertices[2] = vec3<float>(-t, -1.0, 0.0) / denom;
    vertices[3] = vec3<float>(t, -1.0, 0.0) / denom;
    vertices[4] = vec3<float>(1.0, 0.0, t) / denom;
    vertices[5] = vec3<float>(1.0, 0.0, -t) / denom;
    vertices[6] = vec3<float>(-1.0, 0.0, -t) / denom;
    vertices[7] = vec3<float>(-1.0, 0.0, t) / denom;
    vertices[8] = vec3<float>(0.0, t, 1.0) / denom;
    vertices[9] = vec3<float>(0.0, -t, 1.0) / denom;
    vertices[10] = vec3<float>(0.0, -t, -1.0) / denom;
    vertices[11] = vec3<float>(0.0, t, -1.0) / denom;
    triangles.resize(20);
    triangles[0] = vec3<uint32_t>(4, 8, 7);
    triangles[1] = vec3<uint32_t>(4, 7, 9);
    triangles[2] = vec3<uint32_t>(5, 6, 11);
    triangles[3] = vec3<uint32_t>(5, 10, 6);
    triangles[4] = vec3<uint32_t>(0, 4, 3);
    triangles[5] = vec3<uint32_t>(0, 3, 5);
    triangles[6] = vec3<uint32_t>(2, 7, 1);
    triangles[7] = vec3<uint32_t>(2, 1, 6);
    triangles[8] = vec3<uint32_t>(8, 0, 11);
    triangles[9] = vec3<uint32_t>(8, 11, 1);
    triangles[10] = vec3<uint32_t>(9, 10, 3);
    triangles[11] = vec3<uint32_t>(9, 2, 10);
    triangles[12] = vec3<uint32_t>(8, 4, 0);
    triangles[13] = vec3<uint32_t>(11, 0, 5);
    triangles[14] = vec3<uint32_t>(4, 9, 3);
    triangles[15] = vec3<uint32_t>(5, 3, 10);
    triangles[16] = vec3<uint32_t>(7, 8, 1);
    triangles[17] = vec3<uint32_t>(6, 1, 11);
    triangles[18] = vec3<uint32_t>(7, 2, 9);
    triangles[19] = vec3<uint32_t>(6, 10, 2);
    }

  std::vector<jtk::vec3<float>> vertices;
  std::vector<jtk::vec3<uint32_t>> triangles;
  };