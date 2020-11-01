#include "sp.h"
#include "shape_predictor.h"

#include <iostream>

bool read_from_file(sp& shape_pred, const std::string& filename)
  {
  shape_pred.flip_horizontal = false;
  shape_pred.odl = sp::odl_none;
  try
    {
    //shape_pred.p_shape_predictor.reset(new shape_predictor(filename));
    shape_pred.p_shape_predictor = std::make_shared< shape_predictor>(filename);
    }
  catch (...)
    {
    shape_pred.p_shape_predictor.reset();
    return false;
    }
  return true;
  }

std::vector<std::pair<long, long>> predict(const sp& shape_pred, const rect& r, int w, int h, int stride, const uint32_t* p_image)
  {
  std::vector<std::pair<long, long>> landmarks;
  if (shape_pred.p_shape_predictor.get())
    {
    landmarks = shape_pred.p_shape_predictor->predict(r, w, h, stride, p_image, shape_pred.flip_horizontal);
    }
  return landmarks;
  }

void swap(sp& left, sp& right)
  {
  std::swap(left.flip_horizontal, right.flip_horizontal);
  std::swap(left.p_shape_predictor, right.p_shape_predictor);
  std::swap(left.odl, right.odl);
  }

namespace
  {
  std::string link_to_string(sp::object_detector_link lnk)
    {
    switch (lnk)
      {
      case sp::odl_none: return "none";
      case sp::odl_facial: return "facial detector";
      case sp::odl_ear_left: return "ear left detector";
      case sp::odl_ear_right: return "ear right detector";
      }
    return "unknown";
    }
  }

void info(const sp& shape_pred)
  {
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "SHAPE PREDICTOR" << std::endl;
  std::cout << "Flip horizontal: " << (shape_pred.flip_horizontal ? "Yes" : "No") << std::endl;
  std::cout << "Linked: " << link_to_string(shape_pred.odl) << std::endl;
  std::cout << "---------------------------------------" << std::endl;
  }