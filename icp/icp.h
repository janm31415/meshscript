#pragma once

#include <jtk/mat.h>
#include <jtk/vec.h>
#include <vector>

#include <jtk/point_tree.h>


namespace jtk
  {

  class icp
    {
    protected:
      struct point
        {
        vec3<float> pt;
        float& operator [] (size_t i)
          {
          return pt[i];
          }
        float operator [] (size_t i) const
          {
          return pt[i];
          }
        uint64_t idx;
        };

      struct point_tree_traits
        {
        typedef float value_type;
        enum { dimension = 3 };
        typedef point point_type;
        };

    public:
      icp(const std::vector<vec3<float>>& model_points);

      virtual ~icp();

      float fit(matf16& transformation, const std::vector<vec3<float>>& template_points, float inlier_distance = -1, void* data = nullptr);

      float inlier_ratio() const
        {
        return _inlier_ratio;
        }

      uint32_t inlier_count() const
        {
        return (uint32_t)_active.size();
        }

      void set_maximum_iterations(uint32_t max_iter)
        {
        _max_iter = max_iter;
        }

      void set_tolerance(float tolerance)
        {
        _tolerance = tolerance;
        }

      void set_inlier_shrink_factor(float f)
        {
        _inlier_shrink_factor = f;
        }

      void set_minimum_inlier_distance(float d)
        {
        _minimum_inlier_distance = d;
        }

    private:
      virtual std::vector<uint64_t> inliers(const std::vector<vec3<float>>& template_points, const matf16& transformation, float inlier_distance, void* data = nullptr) = 0;
      virtual float residual(const std::vector<vec3<float>>& template_points, const matf16& transformation, const std::vector<uint64_t>& active) = 0;
      virtual float fit_step(matf16& transformation, const std::vector<vec3<float>>& template_points, const std::vector<uint64_t>& active) = 0;

    protected:
      std::vector<point> _model_points;
      point_tree<point_tree_traits> _tree;
      std::vector<uint64_t> _active;
      float _inlier_ratio;
      uint32_t _max_iter;
      float _tolerance;
      float _inlier_shrink_factor;
      float _minimum_inlier_distance;
    };

  }