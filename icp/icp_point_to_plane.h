#pragma once

#include "icp.h"

namespace jtk
  {

  class icp_point_to_plane : public icp
    {
    public:
      icp_point_to_plane(const std::vector<vec3<float>>& model_points, int number_of_neighbours = 10) : icp(model_points)
        {
        _normals = compute_normals(number_of_neighbours);
        }
      virtual ~icp_point_to_plane() {}

    private:
      virtual std::vector<uint64_t> inliers(const std::vector<vec3<float>>& template_points, const matf16& transformation, float inlier_distance, void* data);
      virtual float residual(const std::vector<vec3<float>>& template_points, const matf16& transformation, const std::vector<uint64_t>& active);
      virtual float fit_step(matf16& transformation, const std::vector<vec3<float>>& template_points, const std::vector<uint64_t>& active);

      std::vector<vec3<float>> compute_normals(int number_of_neighbours);

    private:
      std::vector<vec3<float>> _normals;
    };


  }