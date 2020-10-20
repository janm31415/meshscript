#pragma once

#include <memory>
#include <string>
#include <vector>

#include "objects.h"

class shape_predictor;

struct sp
  {
  std::unique_ptr<shape_predictor> p_shape_predictor;
  };

bool read_from_file(sp& shape_pred, const std::string& filename);


std::vector<std::pair<long, long>> predict(const sp& shape_pred, const rect& r, int w, int h, int stride, const uint32_t* p_image);