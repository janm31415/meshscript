#pragma once

#include <jtk/vec.h>
#include <stdint.h>
#include "libpoisson_api.h"

#include <vector>

struct poisson_reconstruction_screened_parameters
  {
  poisson_reconstruction_screened_parameters()
    {
    depth = 8;
    octree_depth = 5;
    samples_per_node = 1.5;
    solver_divide = 0;
    }

  int depth;
  int octree_depth;
  double samples_per_node;
  int solver_divide;
  };

LIBPOISSON_API void poisson_reconstruction_screened(std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& pts, const std::vector<jtk::vec3<float>>& normals, const poisson_reconstruction_screened_parameters& par);

LIBPOISSON_API void poisson_reconstruction_screened(std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertex_colors, const std::vector<jtk::vec3<float>>& pts, const std::vector<jtk::vec3<float>>& normals, const std::vector<uint32_t>& colors, const poisson_reconstruction_screened_parameters& par);
