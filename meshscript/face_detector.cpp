#include "face_detector.h"

#include <dlib/image_processing/frontal_face_detector.h>
#include <dlib/image_processing.h>
#include <dlib/image_io.h>

#include <dlib/image_saver/save_png.h>

#include <vector>
#include <iostream>
#include <streambuf>

#include <jtk/image.h>

struct face_detector_data
  {
  dlib::frontal_face_detector detector;
  };

face_detector::face_detector() : mp_data(new face_detector_data)
  {
  mp_data->detector = dlib::get_frontal_face_detector();
  }

face_detector::~face_detector()
  {
  }

std::vector<rect> face_detector::detect(int w, int h, int stride, const uint32_t* p_image)
  {
  std::vector<rect> out;
  using namespace dlib;
  array2d<rgb_pixel> img;
  img.set_size(h, w);
  for (size_t y = 0; y < h; ++y)
    {
    const uint32_t* p_im = p_image + y * stride;
    for (size_t x = 0; x < w; ++x, ++p_im)
      {
      uint32_t col = *p_im;
      auto& pix = img[y][x];
      pix.red = col & 0xff;
      pix.green = (col >> 8) & 0xff;
      pix.blue = (col >> 16) & 0xff;
      }
    }
  std::vector<rectangle> dets = mp_data->detector(img);
  for (unsigned long j = 0; j < dets.size(); ++j)
    {
    rect r;
    r.x = dets[j].left();
    r.y = dets[j].top();
    r.w = dets[j].width();
    r.h = dets[j].height();
    out.push_back(r);
    }
  return out;
  }

void face_detector::draw_prediction_rgba(int w, int h, int stride, uint32_t* p_image, const std::vector<rect>& rectangles) const
  {
  uint32_t line_color = 0xff00ff00;
  uint32_t* p_im = (p_image);
  for (const auto& r : rectangles)
    {
    jtk::draw_box(p_image, r.x, r.y, r.w, r.h, w, h, stride, line_color);
    }
  }