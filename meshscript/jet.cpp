#include "jet.h"

namespace
  {
  float interpolate(float val, float y0, float x0, float y1, float x1)
    {
    return (val - x0)*(y1 - y0) / (x1 - x0) + y0;
    }

  float base(float val)
    {
    if (val <= -0.75f)
      return 0.f;
    else if (val <= -0.25f)
      return interpolate(val, 0.f, -0.75f, 1.f, -0.25f);
    else if (val <= 0.25f)
      return 1.f;
    else if (val <= 0.75f)
      return interpolate(val, 1.f, 0.25f, 0.f, 0.75f);
    else return 0.f;
    }

  }

float jet_red(float gray) 
  {
  return base(gray - 0.5f);
  }

float jet_green(float gray) 
  {
  return base(gray);
  }

float jet_blue(float gray) 
  {
  return base(gray + 0.5f);
  }

void jet_cold_to_hot(float& r, float& g, float& b, float v)
  {
  const float vmin = 0.f;
  const float vmax = 1.f;
  if (v < vmin)
    v = vmin;
  if (v > vmax)
    v = vmax;
  const float dv = vmax - vmin;

  r = g = b = 1.f;
  if (v < (vmin + 0.25f * dv)) 
    {
    r = 0.f;
    g = 4.f * (v - vmin) / dv;
    }
  else if (v < (vmin + 0.5f * dv)) 
    {
    r = 0.f;
    b = 1.f + 4.f * (vmin + 0.25f * dv - v) / dv;
    }
  else if (v < (vmin + 0.75f * dv)) 
    {
    r = 4.f * (v - vmin - 0.5f * dv) / dv;
    b = 0.f;
    }
  else {
    g = 1.f + 4.f * (vmin + 0.75f * dv - v) / dv;
    b = 0.f;
    }
  }