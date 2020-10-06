#include "icp_point_to_point.h"
#include <jtk/fitting.h>

namespace jtk
  {

  std::vector<uint64_t> icp_point_to_point::inliers(const std::vector<vec3<float>>& template_points, const matf16& transformation, float inlier_distance, void* data)
    {
    std::vector<uint64_t> inliers;
    for (uint64_t i = 0; i < (uint64_t)template_points.size(); ++i)
      {
      matf4 v(4);
      v << template_points[i][0], template_points[i][1], template_points[i][2], 1;
      matf4 v_transf = transformation * v;
      point vt;
      vt.pt = vec3<float>(v_transf(0), v_transf(1), v_transf(2));
      float dist;
      _tree.find_nearest(dist, vt);
      if (dist < inlier_distance)
        inliers.push_back(i);
      }
    return inliers;
    }

  float icp_point_to_point::residual(const std::vector<vec3<float>>& template_points, const matf16& transformation, const std::vector<uint64_t>& active)
    {
    float residual = 0.0;
    for (uint64_t i = 0; i < (uint64_t)active.size(); ++i)
      {
      matf4 v(4);
      v << template_points[active[i]][0], template_points[active[i]][1], template_points[active[i]][2], 1;
      matf4 v_transf = transformation * v;
      point vt;
      vt.pt = vec3<float>(v_transf(0), v_transf(1), v_transf(2));
      float dist;
      _tree.find_nearest(dist, vt);
      residual += dist;
      }
    residual /= (float)active.size();
    return residual;
    }

  float icp_point_to_point::fit_step(matf16& transformation, const std::vector<vec3<float>>& template_points, const std::vector<uint64_t>& active)
    {
    matf source(active.size(), 3);
    matf destination(active.size(), 3);
    for (uint64_t i = 0; i < active.size(); ++i)
      {
      matf4 v(4);
      v << template_points[active[i]][0], template_points[active[i]][1], template_points[active[i]][2], 1;
      matf4 v_transf = transformation * v;
      point vt;
      vt.pt = vec3<float>(v_transf(0), v_transf(1), v_transf(2));
      float dist;
      auto nearest_source = _tree.find_nearest(dist, vt);
      source(i, 0) = nearest_source.pt[0];
      source(i, 1) = nearest_source.pt[1];
      source(i, 2) = nearest_source.pt[2];
      destination(i, 0) = v_transf(0);
      destination(i, 1) = v_transf(1);
      destination(i, 2) = v_transf(2);
      }
    matf16 delta_transformation = npoint(destination, source);

    transformation = delta_transformation * transformation;

    float l2_norm_rotation = norm(block(delta_transformation, 0, 0, 3, 3) - identity(3, 3));
    float l2_norm_translation = norm(block(delta_transformation, 0, 3, 3, 1));
    return std::max<float>(l2_norm_rotation, l2_norm_translation);
    }


  }