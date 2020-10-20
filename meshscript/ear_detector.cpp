#include "ear_detector.h"

#include <dlib/image_processing.h>
#include <dlib/image_io.h>

#include <dlib/image_saver/save_png.h>

#include <vector>
#include <iostream>
#include <streambuf>

#include <jtk/image.h>
#include <jtk/file_utils.h>

struct ear_detector_data
  {
  typedef dlib::scan_fhog_pyramid<dlib::pyramid_down<6> > image_scanner_type;
  dlib::object_detector<image_scanner_type> detector;
  };

ear_detector::ear_detector() : mp_data(new ear_detector_data)
  {
  std::string path = jtk::get_folder(jtk::get_executable_path()) + std::string("data/ear_detector.svm");
  dlib::deserialize(path) >> mp_data->detector;
  }

ear_detector::~ear_detector()
  {
  }

std::vector<rect> ear_detector::detect(int w, int h, int stride, const uint32_t* p_image, bool right_ear)
  {
  std::vector<rect> out;
  using namespace dlib;
  array2d<rgb_pixel> img;
  img.set_size(h, w);
  if (right_ear)
    {
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
    }
  else
    {
    for (size_t y = 0; y < h; ++y)
      {
      const uint32_t* p_im = p_image + y * stride;
      for (size_t x = 0; x < w; ++x, ++p_im)
        {
        uint32_t col = *p_im;
        auto& pix = img[y][w-x-1];
        pix.red = col & 0xff;
        pix.green = (col >> 8) & 0xff;
        pix.blue = (col >> 16) & 0xff;
        }
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
    if (!right_ear)
      {
      r.x = w - r.x - r.w;
      }
    out.push_back(r);
    }
  return out;
  }

void ear_detector::draw_prediction_rgba(int w, int h, int stride, uint32_t* p_image, const std::vector<rect>& rectangles) const
  {
  uint32_t line_color = 0xff00ff00;
  for (const auto& r : rectangles)
    {
    jtk::draw_box(p_image, r.x, r.y, r.w, r.h, w, h, stride, line_color);
    }
  }