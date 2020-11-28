#pragma once

#include <jtk/vec.h>
#include <vector>

struct fill_hole_minimal_surface_parameters
  {
  enum class e_edge_swap_strategy
    {
    minimal_area,
    circum_center
    };

  fill_hole_minimal_surface_parameters();

  int iterations;
  int number_of_rings;
  e_edge_swap_strategy edge_swap_strategy;
  };

void fill_hole_minimal_surface(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<uint32_t>& hole, const fill_hole_minimal_surface_parameters& params);