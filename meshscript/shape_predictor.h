#pragma once

#include <string>

#include <stdint.h>

#include <memory>
#include <vector>

#include "objects.h"

struct shape_predictor_data;

class shape_predictor
  {
  public:
    shape_predictor(const std::string& shape_predictor_filename);
    ~shape_predictor();

    std::vector<std::pair<long, long>> predict(const rect& r, int w, int h, int stride, const uint32_t* p_image);

    void draw_prediction_rgba(int w, int h, int stride, uint32_t* p_image, const std::vector<std::pair<long, long>>& points) const;

  private:
    std::unique_ptr<shape_predictor_data> mp_data;
  };