#pragma once

#include "jtk/jtk/vec.h"
#include <vector>
#include <stdint.h>

std::vector<std::pair<float, float>> lscm(const jtk::vec3<uint32_t>* triangles, uint32_t nr_of_triangles, const jtk::vec3<float>* vertices, uint32_t nr_of_vertices);

void scale_to_unit(std::vector<std::pair<float, float>>& uv);