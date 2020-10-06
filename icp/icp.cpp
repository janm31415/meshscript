#include "icp.h"
#include <numeric>

namespace jtk
  {

  icp::icp(const std::vector<vec3<float>>& model_points) : _max_iter(200),
    _tolerance(1e-4), _inlier_shrink_factor(0.9), _minimum_inlier_distance(0.05)
    {
    _model_points.reserve(model_points.size());
    for (uint64_t i = 0; i < model_points.size(); ++i)
      {
      point p;
      p.pt = model_points[i];
      p.idx = i;
      _model_points.push_back(p);
      }
    _tree.efficient_build_tree(_model_points.begin(), _model_points.end());
    }

  icp::~icp()
    {

    }

  float icp::fit(matf16& transformation, const std::vector<vec3<float>>& template_points, float inlier_distance, void* data)
    {
    assert(transformation.cols() == 4 && transformation.rows() == 4);
    if (_tree.empty())
      throw std::runtime_error("icp: empty kd tree");
    if (template_points.size() < 5)
      throw std::runtime_error("icp: I need at least 5 template points");
    if (inlier_distance <= 0.0)
      {
      _active.resize(template_points.size());
      std::iota(_active.begin(), _active.end(), 0);
      _inlier_ratio = 1.0;
      }
    float delta = 1000.0;
    for (uint32_t iter = 0; iter < _max_iter && delta > _tolerance; ++iter)
      {
      if (inlier_distance > 0.0)
        {
        inlier_distance = std::max<float>(inlier_distance*_inlier_shrink_factor, _minimum_inlier_distance);
        _active = inliers(template_points, transformation, inlier_distance, data);
        _inlier_ratio = (float)_active.size() / template_points.size();
        }
      delta = fit_step(transformation, template_points, _active);
      }
    return residual(template_points, transformation, _active);
    }

  }