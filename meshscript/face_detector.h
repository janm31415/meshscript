#pragma once

#include <string>

#include <stdint.h>

#include <memory>
#include <vector>

struct face_detector_data;

class face_detector
  {
  public:
    face_detector(const std::string& shape_predictor);
    ~face_detector();    

    std::vector<std::pair<long, long>> predict(int w, int h, int stride, const uint32_t* p_image);

    void draw_prediction_rgba(int w, int h, int stride, uint32_t* p_image, const std::vector<std::pair<long, long>>& points) const;    

  private:
    std::unique_ptr<face_detector_data> mp_data;
  };