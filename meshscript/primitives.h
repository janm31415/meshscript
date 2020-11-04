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


struct cube
  {
  cube()
    {
    using namespace jtk;
    vertices.reserve(8);
    triangles.reserve(12);
    vertices.emplace_back((float)0, (float)0, (float)1);
    vertices.emplace_back((float)1, (float)0, (float)1);
    vertices.emplace_back((float)1, (float)1, (float)1);
    vertices.emplace_back((float)0, (float)1, (float)1);
    vertices.emplace_back((float)0, (float)0, (float)0);
    vertices.emplace_back((float)1, (float)0, (float)0);
    vertices.emplace_back((float)1, (float)1, (float)0);
    vertices.emplace_back((float)0, (float)1, (float)0);

    triangles.emplace_back(0, 1, 2);
    triangles.emplace_back(0, 2, 3);
    triangles.emplace_back(7, 6, 5);
    triangles.emplace_back(7, 5, 4);
    triangles.emplace_back(1, 0, 4);
    triangles.emplace_back(1, 4, 5);
    triangles.emplace_back(2, 1, 5);
    triangles.emplace_back(2, 5, 6);
    triangles.emplace_back(3, 2, 6);
    triangles.emplace_back(3, 6, 7);
    triangles.emplace_back(0, 3, 7);
    triangles.emplace_back(0, 7, 4);
    }

  std::vector<jtk::vec3<float>> vertices;
  std::vector<jtk::vec3<uint32_t>> triangles;
  };


struct cylinder
  {
  cylinder()
    {
    double twopi = 2.0* 3.14159265358979323846264338327950288419716939937510;
    using namespace jtk;
    int steps = 128;
    vertices.reserve(steps * 2 + 2);
    triangles.reserve(4 * steps);
    for (int i = 0; i < steps; ++i)
      {
      float cs = std::cos((double)i/(double)steps*twopi);
      float sn = std::sin((double)i / (double)steps*twopi);
      vertices.emplace_back(cs, sn, 0.f);
      vertices.emplace_back(cs, sn, 1.f);      
      }
    vertices.emplace_back(0.f, 0.f, 0.f);
    vertices.emplace_back(0.f, 0.f, 1.f);

    for (int i = 0; i < steps-1; ++i)
      {
      triangles.emplace_back(2 * i, 2 * i + 2, 2*i+1);
      triangles.emplace_back(2 * i + 1, 2 * i + 2, 2*i + 3);

      triangles.emplace_back(2 * steps, 2 * i + 2, 2 * i);
      triangles.emplace_back(2 * steps+1, 2 * i + 1, 2 * i + 3);
      }
    triangles.emplace_back(2 * (steps-1), 0, 2 * (steps - 1) + 1);
    triangles.emplace_back(2 * (steps - 1) + 1, 0, 1);

    triangles.emplace_back(2*steps, 0, 2 * (steps - 1));
    triangles.emplace_back(2 * steps + 1, 2 * (steps - 1) + 1, 1);
    }
  std::vector<jtk::vec3<float>> vertices;
  std::vector<jtk::vec3<uint32_t>> triangles;
  };