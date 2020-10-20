#include "sp.h"
#include "shape_predictor.h"

bool read_from_file(sp& shape_pred, const std::string& filename)
  {
  try
    {
    shape_pred.p_shape_predictor.reset(new shape_predictor(filename));
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
    landmarks = shape_pred.p_shape_predictor->predict(r, w, h, stride, p_image);
    }
  return landmarks;
  }