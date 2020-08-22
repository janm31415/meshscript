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
  dlib::shape_predictor sp;
  };

face_detector::face_detector(const std::string& shape_predictor_filename) : mp_data(new face_detector_data)
  {
  mp_data->detector = dlib::get_frontal_face_detector();
  dlib::deserialize(shape_predictor_filename) >> mp_data->sp;
  }

face_detector::~face_detector()
  {

  }

std::vector<std::pair<long, long>> face_detector::predict(int w, int h, int stride, const uint32_t* p_image)
  {
  std::vector<std::pair<long, long>> out;
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
    full_object_detection shape = mp_data->sp(img, dets[j]);
    for (unsigned i = 0; i < 68; ++i)
      {
      out.push_back(std::pair<long, long>(shape.part(i).x(), shape.part(i).y()));
      }
    }
  return out;
  }

void face_detector::draw_prediction_rgba(int w, int h, int stride, uint32_t* p_image, const std::vector<std::pair<long, long>>& points) const
  {
  uint32_t line_color = 0xff000000 | (uint32_t(0) << 16) | (uint32_t(255) << 8) | uint32_t(0);
  uint32_t* p_im = (p_image);
  size_t nr_of_detections = points.size() / 68;
  for (size_t j = 0; j < nr_of_detections; ++j)
    {
    // Around Chin. Ear to Ear
    for (unsigned long i = 1; i <= 16; ++i)
      {
      auto p0 = points[j * 68 + i];
      auto p1 = points[j * 68 + i - 1];
      jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
      }

    // Line on top of nose
    for (unsigned long i = 28; i <= 30; ++i)
      {
      auto p0 = points[j * 68 + i];
      auto p1 = points[j * 68 + i - 1];
      jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
      }

    // left eyebrow
    for (unsigned long i = 18; i <= 21; ++i)
      {
      auto p0 = points[j * 68 + i];
      auto p1 = points[j * 68 + i - 1];
      jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
      }
    // Right eyebrow
    for (unsigned long i = 23; i <= 26; ++i)
      {
      auto p0 = points[j * 68 + i];
      auto p1 = points[j * 68 + i - 1];
      jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
      }
    // Bottom part of the nose
    for (unsigned long i = 31; i <= 35; ++i)
      {
      auto p0 = points[j * 68 + i];
      auto p1 = points[j * 68 + i - 1];
      jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
      }
    // Line from the nose to the bottom part above
    {
    auto p0 = points[j * 68 + 30];
    auto p1 = points[j * 68 + 35];
    jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
    }

    // Left eye
    for (unsigned long i = 37; i <= 41; ++i)
      {
      auto p0 = points[j * 68 + i];
      auto p1 = points[j * 68 + i - 1];
      jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
      }
    {
    auto p0 = points[j * 68 + 36];
    auto p1 = points[j * 68 + 41];
    jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
    }

    // Right eye
    for (unsigned long i = 43; i <= 47; ++i)
      {
      auto p0 = points[j * 68 + i];
      auto p1 = points[j * 68 + i - 1];
      jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
      }
    {
    auto p0 = points[j * 68 + 42];
    auto p1 = points[j * 68 + 47];
    jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
    }

    // Lips outer part
    for (unsigned long i = 49; i <= 59; ++i)
      {
      auto p0 = points[j * 68 + i];
      auto p1 = points[j * 68 + i - 1];
      jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
      }
    {
    auto p0 = points[j * 68 + 48];
    auto p1 = points[j * 68 + 59];
    jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
    }

    // Lips inside part
    for (unsigned long i = 61; i <= 67; ++i)
      {
      auto p0 = points[j * 68 + i];
      auto p1 = points[j * 68 + i - 1];
      jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
      }
    {
    auto p0 = points[j * 68 + 60];
    auto p1 = points[j * 68 + 67];
    jtk::draw_line(p_im, p0.first, p0.second, p1.first, p1.second, w, h, stride, line_color);
    }
    }
  }