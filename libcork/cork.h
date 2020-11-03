#pragma once

#include "libcork_api.h"

#include <ostream>
#include <vector>

#include <jtk/vec.h>


struct cork_options
  {
  cork_options() : use_parallel(true), debug_folder(nullptr), resolve_all_intersections(false), p_str(nullptr) {}
  bool use_parallel;
  const char* debug_folder;
  bool resolve_all_intersections;
  std::ostream* p_str;
  };

LIBCORK_API bool has_intersections(const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& vertices, const cork_options& options);

LIBCORK_API bool is_closed(const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& vertices);

LIBCORK_API bool is_solid(const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& vertices, const cork_options& options);

LIBCORK_API void compute_union(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles_left, const std::vector<jtk::vec3<float>>& vertices_left, const std::vector<jtk::vec3<uint32_t>>& triangles_right, const std::vector<jtk::vec3<float>>& vertices_right, const cork_options& options);

LIBCORK_API void compute_difference(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles_left, const std::vector<jtk::vec3<float>>& vertices_left, const std::vector<jtk::vec3<uint32_t>>& triangles_right, const std::vector<jtk::vec3<float>>& vertices_right, const cork_options& options);

LIBCORK_API void compute_intersection(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles_left, const std::vector<jtk::vec3<float>>& vertices_left, const std::vector<jtk::vec3<uint32_t>>& triangles_right, const std::vector<jtk::vec3<float>>& vertices_right, const cork_options& options);

LIBCORK_API void compute_symmetric_difference(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles_left, const std::vector<jtk::vec3<float>>& vertices_left, const std::vector<jtk::vec3<uint32_t>>& triangles_right, const std::vector<jtk::vec3<float>>& vertices_right, const cork_options& options);

LIBCORK_API void resolve_intersections(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles_in, const std::vector<jtk::vec3<float>>& vertices_in, const cork_options& options);


struct diagnostics
  {
  uint32_t number_of_vertices;
  uint32_t number_of_triangles;
  uint32_t intersections;
  uint32_t degenerate;
  };

LIBCORK_API diagnostics diagnose(const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& vertices, const cork_options& options);


