#pragma once

#include <vector>
#include <stdint.h>
#include <jtk/vec.h>

bool read_ply(const char* filename, std::vector<jtk::vec3<float>>& pts, std::vector<jtk::vec3<float>>& normals, std::vector<uint32_t>& clrs);
  