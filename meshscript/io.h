#pragma once

#include <vector>
#include <stdint.h>
#include <jtk/vec.h>

bool read_ply(const char* filename, std::vector<jtk::vec3<float>>& pts, std::vector<jtk::vec3<float>>& normals, std::vector<uint32_t>& clrs);


bool read_ply(const char* filename, std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<uint32_t>& clrs);
  