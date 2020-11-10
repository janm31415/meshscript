#include "icp_point_to_plane.h"
#include <jtk/fitting.h>

namespace jtk
  {


  std::vector<uint64_t> icp_point_to_plane::inliers(const std::vector<vec3<float>>& template_points, const matf16& transformation, float inlier_distance, void* data)
    {
    (void*)data;
    std::vector<uint64_t> inliers;
    for (uint64_t i = 0; i < (uint64_t)template_points.size(); ++i)
      {
      matf4 v(4);
      v << template_points[i][0], template_points[i][1], template_points[i][2], 1;
      matf4 v_transf = transformation * v;
      point vt;
      vt.pt = vec3<float>(v_transf(0), v_transf(1), v_transf(2));
      float dist;
      auto p = _tree.find_nearest(dist, vt);

      float dx = p.pt[0];
      float dy = p.pt[1];
      float dz = p.pt[2];
      float nx = _normals[p.idx][0];
      float ny = _normals[p.idx][1];
      float nz = _normals[p.idx][2];

      if (std::abs((vt.pt[0] - dx)*nx + (vt.pt[1] - dy)*ny + (vt.pt[2] - dz)*nz) < inlier_distance)
        inliers.push_back(i);
      }
    return inliers;
    }

  float icp_point_to_plane::residual(const std::vector<vec3<float>>& template_points, const matf16& transformation, const std::vector<uint64_t>& active)
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
      auto p = _tree.find_nearest(dist, vt);
      float dx = p.pt[0];
      float dy = p.pt[1];
      float dz = p.pt[2];
      float nx = _normals[p.idx][0];
      float ny = _normals[p.idx][1];
      float nz = _normals[p.idx][2];

      residual += std::abs((vt.pt[0] - dx)*nx + (vt.pt[1] - dy)*ny + (vt.pt[2] - dz)*nz);
      }
    residual /= (float)active.size();
    return residual;
    }

  float icp_point_to_plane::fit_step(matf16& transformation, const std::vector<vec3<float>>& template_points, const std::vector<uint64_t>& active)
    {
    matf A((uint32_t)active.size(), 6);
    matf b((uint32_t)active.size(), 1);
    for (uint32_t i = 0; i < (uint32_t)active.size(); ++i)
      {
      matf4 v(4);
      v << template_points[active[i]][0], template_points[active[i]][1], template_points[active[i]][2], 1;
      matf4 v_transf = transformation * v;
      point vt;
      vt.pt = vec3<float>(v_transf(0), v_transf(1), v_transf(2));
      float dist;
      auto nearest_source = _tree.find_nearest(dist, vt);
      float dx = nearest_source.pt[0];
      float dy = nearest_source.pt[1];
      float dz = nearest_source.pt[2];
      float nx = _normals[nearest_source.idx][0];
      float ny = _normals[nearest_source.idx][1];
      float nz = _normals[nearest_source.idx][2];

      A(i, 0) = nz * vt.pt[1] - ny * vt.pt[2];
      A(i, 1) = nx * vt.pt[2] - nz * vt.pt[0];
      A(i, 2) = ny * vt.pt[0] - nx * vt.pt[1];
      A(i, 3) = nx;
      A(i, 4) = ny;
      A(i, 5) = nz;
      b(i, 0) = nx * dx + ny * dy + nz * dz - nx * vt.pt[0] - ny * vt.pt[1] - nz * vt.pt[2];
      }

    matf36 A_ = transpose(A)*A;
    matf36 b_ = transpose(A)*b;

    solve(A_, b_);

    matf9 R_ = identity(3, 3);
    R_(0, 1) = -b_(2, 0);
    R_(1, 0) = +b_(2, 0);
    R_(0, 2) = +b_(1, 0);
    R_(2, 0) = -b_(1, 0);
    R_(1, 2) = -b_(0, 0);
    R_(2, 1) = +b_(0, 0);

    matf9 W, V;
    svd(R_, W, V);

    matf9 R = R_ * transpose(V);

    if (determinant(R) < 0)
      {
      matf9 B = identity(3, 3);
      B(2, 2) = determinant(R);
      R = V * B * transpose(R_);
      }

    matf16 delta_transformation(4, 4);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        delta_transformation(i, j) = R(i, j);
    delta_transformation(0, 3) = b_(3, 0);
    delta_transformation(1, 3) = b_(4, 0);
    delta_transformation(2, 3) = b_(5, 0);
    delta_transformation(3, 0) = 0;
    delta_transformation(3, 1) = 0;
    delta_transformation(3, 2) = 0;
    delta_transformation(3, 3) = 1;


    transformation = delta_transformation * transformation;

    float l2_norm_rotation = norm(block(delta_transformation, 0, 0, 3, 3) - identity(3, 3));
    float l2_norm_translation = norm(block(delta_transformation, 0, 3, 3, 1));
    return std::max<float>(l2_norm_rotation, l2_norm_translation);
    }

  std::vector<vec3<float>> icp_point_to_plane::compute_normals(int number_of_neighbours)
    {
    std::vector<vec3<float>> normals(_model_points.size());
    for (uint64_t i = 0; i < _model_points.size(); ++i)
      {
      std::vector<point> neighbours = _tree.find_k_nearest(number_of_neighbours, _model_points[i]);
      matf P(neighbours.size(), 3);
      matf3 mu = zeros(1, 3);
      for (int j = 0; j < neighbours.size(); ++j)
        {
        P(j, 0) = neighbours[j].pt[0];
        P(j, 1) = neighbours[j].pt[1];
        P(j, 2) = neighbours[j].pt[2];
        mu(0, 0) += neighbours[j].pt[0];
        mu(0, 1) += neighbours[j].pt[1];
        mu(0, 2) += neighbours[j].pt[2];
        }
      mu = mu / float(neighbours.size());
      matf Q = P - ones<float>(neighbours.size(), 1)*mu;
      matf9 H = transpose(Q)*Q;
      matf9 W, V;
      svd(H, W, V);
      vec3<float> normal(H(0, 2), H(1, 2), H(2, 2));
      normals[_model_points[i].idx] = normal;
      }
    return normals;
    }

  }