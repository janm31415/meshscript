#pragma once

#include "icp.h"

namespace jtk
  {

  class icp_point_to_point : public icp
    {
    public:
      icp_point_to_point(const std::vector<vec3<float>>& model_points) : icp(model_points) {}
      virtual ~icp_point_to_point() {}

    private:
      virtual std::vector<uint64_t> inliers(const std::vector<vec3<float>>& template_points, const matf16& transformation, float inlier_distance, void* data);
      virtual float residual(const std::vector<vec3<float>>& template_points, const matf16& transformation, const std::vector<uint64_t>& active);
      virtual float fit_step(matf16& transformation, const std::vector<vec3<float>>& template_points, const std::vector<uint64_t>& active);

    };


  }