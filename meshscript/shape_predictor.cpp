#include "shape_predictor.h"

#include <dlib/image_processing.h>
#include <dlib/image_io.h>

#include <dlib/image_saver/save_png.h>

#include <vector>
#include <iostream>
#include <streambuf>

#include <jtk/image.h>

struct shape_predictor_data
  {
  dlib::shape_predictor sp;
  };

shape_predictor::shape_predictor(const std::string& shape_predictor_filename) : mp_data(new shape_predictor_data)
  {
  dlib::deserialize(shape_predictor_filename) >> mp_data->sp;
  }

shape_predictor::~shape_predictor()
  {
  }

std::vector<std::pair<long, long>> shape_predictor::predict(const rect& r, int w, int h, int stride, const uint32_t* p_image, bool flip_horizontal)
  {
  std::vector<std::pair<long, long>> out;
  using namespace dlib;
  array2d<rgb_pixel> img;
  img.set_size(h, w);
  if (flip_horizontal)
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
  else
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
  rectangle re;
  if (flip_horizontal)    
    re = rectangle(w - r.x - r.w, r.y, w - r.x - 1, r.y + r.h - 1);
  else
    re = rectangle(r.x, r.y, r.x + r.w - 1, r.y + r.h - 1);

  full_object_detection shape = mp_data->sp(img, re);
  for (unsigned i = 0; i < shape.num_parts(); ++i)
    {
    if (flip_horizontal)
      out.push_back(std::pair<long, long>(w - shape.part(i).x() - 1, shape.part(i).y()));
    else
      out.push_back(std::pair<long, long>(shape.part(i).x(), shape.part(i).y()));
    }
    
  return out;
  }

void shape_predictor::draw_prediction_rgba(int w, int h, int stride, uint32_t* p_image, const std::vector<std::pair<long, long>>& points) const
  {
  uint32_t line_color = 0xff00ff00;
  for (const auto& pt : points)
    {
    jtk::draw_box(p_image, pt.first - 1, pt.second - 1, 3, 3, w, h, stride, line_color);
    }
  }