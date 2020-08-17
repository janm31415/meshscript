#pragma once

#include "morphable_model.h"
#include <jtk/qbvh.h>

struct mm
  {
  jtk::morphable_model m;
  jtk::float4x4 cs;
  bool visible;
  std::vector<jtk::vec3<float>> vertices;
  std::vector<float> coefficients;
  };

bool read_from_file(mm& morph, const std::string& filename);

bool vertices_to_csv(const mm& m, const std::string& filename);
bool triangles_to_csv(const mm& m, const std::string& filename);

bool write_to_file(const mm& morph, const std::string& filename);