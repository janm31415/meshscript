#include "fill_holes.h"

#include <jtk/geometry.h>
#include <jtk/fitting.h>

#include <cassert>
#include <iostream>

/*
Constructing Triangular Meshes of Minimal Area
Wenyu Chen1, Yiyu Cai2, and Jianmin Zheng3.
Computer-Aided Design and Applications, 5(1-4), 2008, 508-518.
*/

fill_hole_minimal_surface_parameters::fill_hole_minimal_surface_parameters()
  {
  number_of_rings = 12;
  iterations = 100;
  }

void fill_hole_minimal_surface(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<uint32_t>& hole, const fill_hole_minimal_surface_parameters& params)
  {
  int s = params.number_of_rings;
  if (s < 1)
    s = 1;

  jtk::vec3<float> center(0.f, 0.f, 0.f);
  for (auto i : hole)
    {
    center = center + vertices[i];
    }
  center = center / (float)hole.size();

  uint32_t center_id = (uint32_t)vertices.size();

  vertices.reserve(vertices.size() + 1 + hole.size()*s);

  vertices.push_back(center);

  for (auto i : hole)
    {
    for (int j = 1; j <= s; ++j)
      {
      jtk::vec3<float> Pij = vertices[i] + (float)j / (float)(s + 1)*(center - vertices[i]);
      vertices.push_back(Pij);
      }
    }

  size_t new_triangles = hole.size()* (s * 2 + 1);

  triangles.reserve(triangles.size() + new_triangles);
  for (uint32_t k0 = 0; k0 < (uint32_t)hole.size(); ++k0)
    {
    uint32_t k1 = (k0 + 1) % hole.size();
    uint32_t i0 = hole[k0];
    uint32_t i1 = hole[k1];
    triangles.emplace_back(i0, center_id + k0 * s + 1, center_id + k1 * s + 1);
    triangles.emplace_back(i1, i0, center_id + k1 * s + 1);
    for (int j = 1; j < s; ++j)
      {
      triangles.emplace_back(center_id + k0 * s + j, center_id + k0 * s + j + 1, center_id + k1 * s + j + 1);
      triangles.emplace_back(center_id + k1 * s + j, center_id + k0 * s + j, center_id + k1 * s + j + 1);
      }
    triangles.emplace_back(center_id + k1 * s + s, center_id + k0 * s + s, center_id);
    }



  jtk::mutable_adjacency_list adj_list((uint32_t)vertices.size(), triangles.data(), (uint32_t)triangles.size());
  for (int iter = 0; iter < params.iterations; ++iter)
    {
    //first laplace fairing
    std::vector<jtk::vec3<float>> vert(vertices.begin() + center_id, vertices.end());
    auto get_vertex = [&](uint32_t i)->jtk::vec3<float>
      {
      if (i < center_id)
        return vertices[i];
      else
        return vert[i - center_id];
      };
    for (uint32_t v = center_id; v < (uint32_t)vertices.size(); ++v)
      { 
      auto onering = jtk::one_ring_vertices_from_vertex(v, adj_list, triangles.data());
      float total_area = 0.f;
      jtk::vec3<float> umbrella(0, 0, 0);
      for (auto v2 : onering)
        {
        auto tria = jtk::triangle_indices_from_edge(v, v2, adj_list);
        assert(tria.size() == 2);
        
        float tria1_area = jtk::length(jtk::cross(get_vertex(triangles[tria[0]][1]) - get_vertex(triangles[tria[0]][0]), get_vertex(triangles[tria[0]][2]) - get_vertex(triangles[tria[0]][0])));
        float tria2_area = jtk::length(jtk::cross(get_vertex(triangles[tria[1]][1]) - get_vertex(triangles[tria[1]][0]), get_vertex(triangles[tria[1]][2]) - get_vertex(triangles[tria[1]][0])));

        float weight = tria1_area + tria2_area;
        //weight = 1.f;
        umbrella = umbrella + weight * get_vertex(v2);

        total_area += weight;
        }
      vertices[v] = umbrella / total_area;
      }

    size_t edge_swaps = 0;
    // now edge swapping
    bool edge_swapped = true;
    while (edge_swapped && edge_swaps < new_triangles*3)
      {      
      edge_swapped = false;
      for (uint32_t v = center_id; v < (uint32_t)vertices.size(); ++v)
        {
        auto onering = jtk::one_ring_vertices_from_vertex(v, adj_list, triangles.data());
        for (auto v2 : onering)
          {
          if (v2 < v)
            {
            auto tria = jtk::triangle_indices_from_edge(v, v2, adj_list);
            if (tria.size() != 2)
              continue;

            auto t1 = triangles[tria[0]];
            auto t2 = triangles[tria[1]];
            int t1v = 0;
            while (t1[t1v] == v || t1[t1v] == v2)
              ++t1v;
            int t2v = 0;
            while (t2[t2v] == v || t2[t2v] == v2)
              ++t2v;
            
            float tria1_area = jtk::length(jtk::cross(vertices[v] - vertices[t1[t1v]], vertices[v2] - vertices[t1[t1v]]));
            float tria2_area = jtk::length(jtk::cross(vertices[v] - vertices[t2[t2v]], vertices[v2] - vertices[t2[t2v]]));
            
            float tria1_area_swapped = jtk::length(jtk::cross(vertices[t1[t1v]] - vertices[v],  vertices[t2[t2v]] - vertices[v]));
            float tria2_area_swapped = jtk::length(jtk::cross(vertices[t1[t1v]] - vertices[v2], vertices[t2[t2v]] - vertices[v2]));

            if (tria1_area_swapped + tria2_area_swapped < tria1_area + tria2_area)
              {
              bool success = jtk::edge_swap(v, v2, triangles, adj_list);             
              edge_swaps |= success;
              if (success)
                ++edge_swaps;
              }
            }
          }
        }

      }
    } // for iter
  }