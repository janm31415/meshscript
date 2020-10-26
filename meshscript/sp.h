#pragma once

#include <memory>
#include <string>
#include <vector>

#include "objects.h"

class shape_predictor;

struct sp
  {
  enum object_detector_link
    {
    odl_none,
    odl_facial,
    odl_ear_left,
    odl_ear_right
    };

  std::unique_ptr<shape_predictor> p_shape_predictor;
  bool flip_horizontal;
  object_detector_link odl;
  };

bool read_from_file(sp& shape_pred, const std::string& filename);


std::vector<std::pair<long, long>> predict(const sp& shape_pred, const rect& r, int w, int h, int stride, const uint32_t* p_image);

void swap(sp& left, sp& right);

void info(const sp& shape_pred);