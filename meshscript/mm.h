#pragma once

#include "morphable_model.h"
#include <jtk/qbvh.h>

struct mesh;

struct mm
  {
  jtk::morphable_model shape, color;
  jtk::float4x4 cs;
  bool visible;
  std::vector<jtk::vec3<float>> vertices;
  std::vector<jtk::vec3<float>> vertex_colors;
  std::vector<float> coefficients, color_coefficients;
  };

bool read_from_file(mm& morph, const std::string& filename);

bool vertices_to_csv(const mm& m, const std::string& filename);
bool triangles_to_csv(const mm& m, const std::string& filename);

bool write_to_file(const mm& morph, const std::string& filename);

void clamp_vertex_colors(std::vector<jtk::vec3<float>>& vertex_colors);

void fit_to_mesh(mm& morph, const mesh& m);

void fit_to_partial_positions(mm& morph, const std::vector<uint32_t>& vertex_indices, const std::vector<jtk::vec3<float>>& vertex_positions);

void fit(mm& morph, const mesh& m, const std::vector<uint32_t>& vertex_indices, const std::vector<jtk::vec3<float>>& vertex_positions);
