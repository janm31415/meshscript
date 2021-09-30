#pragma once

#include "jtk/deformation.h"
#include "jtk/qbvh.h"
#include <memory>
#include <variant>
#include <vector>

struct deform_tool
  {
  std::vector<jtk::vec3<float>> vertices;
  std::vector<jtk::vec3<uint32_t>> triangles;
  std::variant<std::unique_ptr<jtk::warping_tool>, std::unique_ptr<jtk::pushpull_tool>> tool;
  jtk::float4x4 cs = jtk::get_identity();
  bool visible = true;
  float decay_factor = 1.f;
  bool signed_distance = true;
  uint32_t discretization = 50;
  };

void info(const deform_tool& dt);