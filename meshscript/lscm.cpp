#include "lscm.h"

#include "jtk/geometry.h"
#include "jtk/mat.h"
#include "jtk/timer.h"

#include <map>
#include <queue>
#include <iostream>

namespace
  {
  void _automatically_detect_correspondences(std::map<uint32_t, std::pair<float, float>>& correspondences, const jtk::adjacency_list& adj_list, const jtk::vec3<uint32_t>* triangles, uint32_t nr_of_triangles, const jtk::vec3<float>* vertices, uint32_t nr_of_vertices)
    {
    uint32_t first_vertex, second_vertex;
    for (uint32_t v = 0; v < nr_of_vertices; ++v)
      {
      if (jtk::is_boundary_vertex(v, adj_list, triangles))
        {
        first_vertex = v;
        break;
        }
      }

    std::vector<float> distances(nr_of_vertices, -1.f);
    std::queue<uint32_t> qu;
    qu.push(first_vertex);

    distances[first_vertex] = 0.f;
    while (!qu.empty())
      {
      auto vertex = qu.front();
      qu.pop();
      double current_distance = distances[vertex];
      auto one_ring = jtk::one_ring_vertices_from_vertex(vertex, adj_list, triangles);
      for (auto v : one_ring)
        {
        if (jtk::is_boundary_vertex(v, adj_list, triangles))
          {
          float tmp_dist = current_distance + distance(vertices[vertex], vertices[v]);
          float old_dist = distances[v];
          if (old_dist == -1.f || tmp_dist < old_dist)
            {
            qu.push(v);
            distances[v] = tmp_dist;
            }
          }
        }
      } // while (!qu.empty())

    float max_distance = 0.f;

    for (size_t i = 0; i < distances.size(); ++i)
      {
      if (distances[i] > max_distance)
        {
        max_distance = distances[i];
        second_vertex = i;
        }
      }

    correspondences.clear();
    correspondences[first_vertex] = std::pair<float, float >(0.f, 0.f);
    correspondences[second_vertex] = std::pair<float, float >(1.f, 1.f);
    }
  } // anonymous namespace

std::vector<std::pair<float, float>> lscm(const jtk::vec3<uint32_t>* triangles, uint32_t nr_of_triangles, const jtk::vec3<float>* vertices, uint32_t nr_of_vertices)
  {
  jtk::adjacency_list adj_list(nr_of_vertices, triangles, nr_of_triangles);

  std::map<uint32_t, std::pair<float, float>> correspondences;
  _automatically_detect_correspondences(correspondences, adj_list, triangles, nr_of_triangles, vertices, nr_of_vertices);

  std::vector<uint32_t> indices(nr_of_vertices, 0xffffffff);

  uint32_t index = 0;
  for (uint32_t i = 0; i < nr_of_vertices; ++i)
    {
    auto it = correspondences.find(i);
    if (it == correspondences.end())
      indices[i] = index++;
    }

  uint32_t num_of_correspondences = (uint32_t)correspondences.size();

  std::vector<std::pair<float, float>> res;

  uint32_t n = nr_of_triangles;
  uint32_t m = nr_of_vertices - num_of_correspondences;

  jtk::smatf A(2 * n, 2 * m);
  jtk::matf b = jtk::zeros<float>(2 * n, 1);

  for (uint32_t i = 0; i < nr_of_triangles; ++i)
    {
    auto tria = triangles[i];
    jtk::vec3<float> points[3] = { vertices[tria[0]], vertices[tria[1]], vertices[tria[2]] };
    auto normal = jtk::cross(points[1] - points[0], points[2] - points[0]);
    float area2 = jtk::length(normal);
    normal = jtk::normalize(normal);

    // Local triangle basis.
    auto e1 = jtk::normalize(points[1] - points[0]);
    auto e2 = jtk::cross(normal, e1);

    // Vertices coordinates in the local triangle basis.
    float x1 = 0.f;
    float y1 = 0.f;
    float x2 = jtk::length(points[1] - points[0]);
    float y2 = 0.f;
    float x3 = jtk::dot(points[2] - points[0], e1);
    float y3 = jtk::dot(points[2] - points[0], e2);

    float w_real[3] = { x3 - x2, x1 - x3, x2 - x1 };
    float w_image[3] = { y3 - y2, y1 - y3, y2 - y1 };
    float sqrt_dt = std::sqrt(area2);

    for (int j = 0; j < 3; ++j)
      {
      if (indices[tria[j]] != 0xffffffff)
        {
        A.put(i, indices[tria[j]]) = w_real[j] / sqrt_dt;
        A.put(i+n, indices[tria[j]] + m) = w_real[j] / sqrt_dt;

        A.put(i + n, indices[tria[j]]) = w_image[j] / sqrt_dt;
        A.put(i, indices[tria[j]] + m) = -w_image[j] / sqrt_dt;
        }
      else
        {
        auto uv = correspondences[tria[j]];
        b(i) -= (w_real[j] / sqrt_dt)*uv.first;
        b(i) -= -(w_image[j] / sqrt_dt)*uv.second;

        b(i + n) -= (w_image[j] / sqrt_dt)*uv.first;
        b(i + n) -= (w_real[j] / sqrt_dt)*uv.second;
        }
      }
    }

  auto At = transpose(A);
  auto AtA = At * A;
  jtk::matf Atb;
  Atb.noalias() = At * b;

  jtk::timer t;

  t.start();
  jtk::diagonal_preconditioner<float> P(AtA);
  std::cout << "Preconditioner calculated! -> " << t.time_elapsed() << "s" << std::endl;
  t.start();
  jtk::symmetric_sparse_matrix_wrapper<float> AtA_lower(AtA);
  std::cout << "Symmetric sparse wrapper calculated! -> " << t.time_elapsed() << "s" << std::endl;

  jtk::matf solution;
  float residu;
  uint32_t iterations;
  jtk::matf x0 = jtk::zeros(AtA.rows(), 1);
  auto tolerance = std::numeric_limits<float>::epsilon();

  t.start();
  jtk::preconditioned_conjugate_gradient(solution, residu, iterations, AtA_lower, P, Atb, x0, tolerance);
  std::cout << "solution calculated! -> " << t.time_elapsed() << "s" << std::endl;
  std::cout << "iterations: " << iterations << std::endl;
  std::cout << "residu: " << residu << std::endl;
  
  for (uint32_t index = 0; index < nr_of_vertices; ++index)
    {
    if (indices[index] != -1)
      {
      std::pair<float, float> uv_coord;
      uv_coord.first = solution(indices[index]);
      uv_coord.second = solution(indices[index] + m);
      res.push_back(uv_coord);
      }
    else
      {
      res.push_back(correspondences[index]);
      }
    }

  return res;
  }