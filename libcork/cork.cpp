#include "cork.h"
#include "ga4.h"
#include "isct.h"
#include "fixisct.h"

#include "quantization.h"

#include "aabvh.h"

#include <string.h>

#include <jtk/concurrency.h>
#include <jtk/geometry.h>
#include <jtk/timer.h>

#include "union_find.h"

#include <map>
#include <queue>
#include <numeric>

#include <stdint.h>

#include <array>
#include <atomic>
#include <algorithm>
#include <mutex>

//#define REAL double

extern "C"
  {
#include "triangle.h"
  }

#define EDGE_N 20

#undef min
#undef max

namespace
  {

  std::vector<jtk::vec3<double>> vertices_to_double(const std::vector<jtk::vec3<float>>& in)
    {
    std::vector<jtk::vec3<double>> out;
    out.reserve(in.size());
    for (const auto& pt : in)
      {
      jtk::vec3<double> p((double)pt[0], (double)pt[1], (double)pt[2]);
      out.push_back(p);
      }
    return out;
    }


  std::vector<jtk::vec3<float>> vertices_to_float(const std::vector<jtk::vec3<double>>& in)
    {
    std::vector<jtk::vec3<float>> out;
    out.reserve(in.size());
    for (const auto& pt : in)
      {
      jtk::vec3<float> p((float)pt[0], (float)pt[1], (float)pt[2]);
      out.push_back(p);
      }
    return out;
    }

  struct triangle_info
    {
    uint8_t bit;
    uint32_t component_id;
    uint32_t parent;
    bool marked;
    };

  template <class _Type, class TFunctor>
  void _parallel_for(_Type first, _Type last, TFunctor fun, const cork_options& options)
    {
    if (options.use_parallel)
      jtk::parallel_for(first, last, fun);
    else
      {
      for (_Type i = first; i != last; ++i)
        fun(i);
      }
    }

  std::string _get_debug_file(const std::string& filename, const cork_options& options)
    {
    if (!options.debug_folder)
      return std::string();
    std::string path = std::string(options.debug_folder) + "\\" + filename;
    return path;
    }

  template <typename T, class Compare = std::less<T>>
  std::vector<uint32_t> sort_indices(const std::vector<T>& v, Compare comp = Compare())
    {
    std::vector<uint32_t> indices((uint32_t)v.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](uint32_t i, uint32_t j) { return comp(v[i], v[j]); });
    return indices;
    }

  double get_magnitude(const std::vector<jtk::vec3<double>>& vertices)
    {
    double max_magnitude = 0.0;
    for (const auto& v : vertices)
      {
      max_magnitude = std::max<double>(max_magnitude, std::abs(v[0]));
      max_magnitude = std::max<double>(max_magnitude, std::abs(v[1]));
      max_magnitude = std::max<double>(max_magnitude, std::abs(v[2]));
      }
    return max_magnitude;
    }

  std::vector<jtk::vec3<double>> get_scaled_vertices(const std::vector<jtk::vec3<double>>& vertices, const quantization& q)
    {
    std::vector<jtk::vec3<double>> quantized(vertices.size());
    uint32_t ind = 0;
    for (auto v : vertices)
      {
      quantized[ind][0] = q(v[0]);
      quantized[ind][1] = q(v[1]);
      quantized[ind][2] = q(v[2]);
      ++ind;
      }
    return quantized;
    }

  inline double drand(double min, double max)
    {
    const double invMAX = 1.0 / double(RAND_MAX);
    double rand0to1 = double(std::rand())*invMAX;
    return (max - min)*rand0to1 + min;
    }

  void perturb_vertices(std::vector<jtk::vec3<double>>& vertices, const quantization& q)
    {
    const double eps = q.get_reshrink_factor()*100.0;
    for (auto iter = vertices.rbegin(); iter != vertices.rend(); ++iter)
      {
      auto& v = *iter;
      jtk::vec3<double> perturb(q(drand(-eps, eps)), q(drand(-eps, eps)), q(drand(-eps, eps)));
      v = v + perturb;
      }
    }

  jtk::boundingbox3d<double> compute_bb(const jtk::vec3<uint32_t>& t, const std::vector<jtk::vec3<double>>& vertices)
    {
    auto bb = make_boundingbox(vertices[t[0]], vertices[t[1]]);
    bb.min = jtk::min(bb.min, vertices[t[2]]);
    bb.max = jtk::max(bb.max, vertices[t[2]]);
    return bb;
    }

  jtk::boundingbox3d<double> compute_bb(const jtk::vec2<uint32_t>& e, const std::vector<jtk::vec3<double>>& vertices)
    {
    auto bb = make_boundingbox(vertices[e[0]], vertices[e[1]]);
    return bb;
    }

  /*
  edges_to_triangles_map

  Edges are stored by their corresponding vertex indices.
  The convention is that for and edge (v0, v1) we always have that v0 < v1.
  If edge_map is of type std::vector<std::map<uint32_t, std::vector<uint32_t> > >  and obtained from this method, then we can
  access the triangle indices of the triangles that have edge (v0, v1) as follows:
  std::vector<uint32_t> triangle_indices = edge_map[v0][v1];
  if edge_map is a const object, then we have to use the somewhat less readable version:
  std::vector<uint32_t> triangle_indices = edge_map[v0].find(v1)->second;
  Note that the convention v0 < v1 is important. If we swap positions of v0 and v1 as below (still with v0 < v1), then the call
  std::vector<uint32_t> triangle_indices = edge_map[v1].find(v0)->second;
  will probably crash, while the non-const call
  std::vector<uint32_t> triangle_indices = edge_map[v1][v0];
  will return an incorrect empty list.
  */
  typedef std::vector<std::map<uint32_t, std::vector<uint32_t> > > e_t_map;
  e_t_map edges_to_triangles_map(uint32_t nr_of_vertices, const std::vector<jtk::vec3<uint32_t>>& triangles)
    {
    typedef std::pair<uint32_t, uint32_t> edge_type;
    e_t_map edge_map(nr_of_vertices);
    for (uint32_t t = 0; t < (uint32_t)triangles.size(); ++t)
      {
      auto tria = triangles[t];
      for (uint32_t i = 0; i < 3; ++i)
        {
        uint32_t j = (i + 1) % 3;
        edge_type e(tria[i] < tria[j] ? tria[i] : tria[j], tria[i] < tria[j] ? tria[j] : tria[i]);
        edge_map[e.first][e.second].push_back(t);
        }
      }
    return edge_map;
    }

  aabvh<jtk::vec2<uint32_t>, EDGE_N> compute_edge_aabvh(const e_t_map& edge_map, const std::vector<jtk::vec3<double>>& vertices)
    {
    std::vector<aabvh_entity<jtk::vec2<uint32_t>>> edges;
    uint32_t ind = 0;
    for (const auto& e1 : edge_map)
      {
      for (const auto& e2 : e1)
        {
        jtk::vec2<uint32_t> ed(ind, e2.first);
        aabvh_entity<jtk::vec2<uint32_t>> ent;
        ent.bbox = compute_bb(ed, vertices);
        ent.point = (ent.bbox.min + ent.bbox.max) / 2.0;
        ent.idx = ed;
        edges.push_back(ent);
        }
      ++ind;
      }
    aabvh<jtk::vec2<uint32_t>, EDGE_N> edge_aabvh(edges);
    return edge_aabvh;
    }

  aabvh<jtk::vec2<uint32_t>, EDGE_N> compute_edge_aabvh(const e_t_map& edge_map, const std::vector<jtk::vec3<double>>& vertices, const std::vector<triangle_info>& info, int mesh_id)
    {
    std::vector<aabvh_entity<jtk::vec2<uint32_t>>> edges;
    uint32_t ind = 0;
    for (const auto& e1 : edge_map)
      {
      for (const auto& e2 : e1)
        {
        if (mesh_id >= 0 && (info[e2.second.front()].bit & 1) != mesh_id) // we consider only edges of first mesh
          break;
        jtk::vec2<uint32_t> ed(ind, e2.first);
        aabvh_entity<jtk::vec2<uint32_t>> ent;
        ent.bbox = compute_bb(ed, vertices);
        ent.point = (ent.bbox.min + ent.bbox.max) / 2.0;
        ent.idx = ed;
        edges.push_back(ent);
        }
      ++ind;
      }
    aabvh<jtk::vec2<uint32_t>, EDGE_N> edge_aabvh(edges);
    return edge_aabvh;
    }

  aabvh<uint32_t> compute_triangle_aabvh(const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<double>>& vertices)
    {
    std::vector<aabvh_entity<uint32_t>> tr;
    for (uint32_t t = 0; t < (uint32_t)triangles.size(); ++t)
      {
      aabvh_entity<uint32_t> ent;
      auto tria = triangles[t];
      ent.bbox = compute_bb(tria, vertices);
      ent.point = (vertices[tria[0]] + vertices[tria[1]] + vertices[tria[2]]) / 3.0;
      ent.idx = t;
      tr.push_back(ent);
      }
    aabvh<uint32_t> triangle_aabvh(tr);
    return triangle_aabvh;
    }

  void aabvh_edge_triangle(std::function<bool(const jtk::vec2<uint32_t>& e, const jtk::vec3<uint32_t>& t, uint32_t tria_index)> func, const aabvh<jtk::vec2<uint32_t>, EDGE_N> & edge_aabvh, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<double>>& vertices)
    {
    bool aborted = false;
    uint32_t tria_index = 0;
    for (const auto& t : triangles)
      {
      if (!aborted)
        {
        auto bbox = compute_bb(t, vertices);
        edge_aabvh.for_each_in_box(bbox, [&](const jtk::vec2<uint32_t>& ed)
          {
          if (!func(ed, t, tria_index))
            aborted = true;
          });
        }
      ++tria_index;
      }
    }

  void aabvh_edge_triangle_parallel(std::function<bool(const jtk::vec2<uint32_t>& e, const jtk::vec3<uint32_t>& t, uint32_t tria_index)> func, const aabvh<jtk::vec2<uint32_t>, EDGE_N> & edge_aabvh, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<double>>& vertices, const cork_options& options)
    {
    std::atomic<bool> aborted{ false };
    _parallel_for((uint32_t)0, (uint32_t)triangles.size(), [&](uint32_t tria_index)
      {
      if (!aborted)
        {
        auto bbox = compute_bb(triangles[tria_index], vertices);
        edge_aabvh.for_each_in_box(bbox, [&](const jtk::vec2<uint32_t>& ed)
          {
          if (!func(ed, triangles[tria_index], tria_index))
            aborted = true;
          });
        }
      }, options);
    }

  void aabvh_edge_triangle_parallel_allow_degeneracies(std::function<bool(const jtk::vec2<uint32_t>& e, const jtk::vec3<uint32_t>& t, uint32_t tria_index)> func, const aabvh<jtk::vec2<uint32_t>, EDGE_N> & edge_aabvh, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<double>>& vertices, const cork_options& options)
    {
    _parallel_for((uint32_t)0, (uint32_t)triangles.size(), [&](uint32_t tria_index)
      {
        auto bbox = compute_bb(triangles[tria_index], vertices);
        edge_aabvh.for_each_in_box(bbox, [&](const jtk::vec2<uint32_t>& ed)
          {
          func(ed, triangles[tria_index], tria_index);
          });        
      }, options);
    }

  void aabvh_edge_triangle_parallel(std::function<bool(const jtk::vec2<uint32_t>& e, const jtk::vec3<uint32_t>& t, uint32_t tria_index)> func, const aabvh<jtk::vec2<uint32_t>, EDGE_N> & edge_aabvh, const std::vector<triangle_info>& info, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<double>>& vertices, int mesh_id, const cork_options& options)
    {
    std::atomic<bool> aborted{ false };
    _parallel_for((uint32_t)0, (uint32_t)triangles.size(), [&](uint32_t tria_index)
      {
      if (!aborted && (mesh_id < 0 || (info[tria_index].bit & 1) != mesh_id))
        {
        auto bbox = compute_bb(triangles[tria_index], vertices);
        edge_aabvh.for_each_in_box(bbox, [&](const jtk::vec2<uint32_t>& ed)
          {
          if (!func(ed, triangles[tria_index], tria_index))
            aborted = true;
          });
        }
      }, options);
    }

  inline bool has_common_vertex(const jtk::vec2<uint32_t>& e, const jtk::vec3<uint32_t>& t)
    {
    return (e[0] == t[0] ||
      e[0] == t[1] ||
      e[0] == t[2] ||
      e[1] == t[0] ||
      e[1] == t[1] ||
      e[1] == t[2]);
    }

  inline int64_t get_common_vertex(const jtk::vec3<uint32_t>& t0, const jtk::vec3<uint32_t>& t1)
    {
    for (uint32_t i = 0; i < 3; ++i)
      {
      for (uint32_t j = 0; j < 3; ++j)
        {
        if (t0[i] == t1[j])
          return int64_t(t0[i]);
        }
      }
    return -1;
    }

  ga::vector<double> make_ga_vector(const jtk::vec3<double>& pt)
    {
    ga::vector<double> v;
    v.e0 = pt[0];
    v.e1 = pt[1];
    v.e2 = pt[2];
    v.e3 = 1.0;
    return v;
    }

  ga::vector<fixint::limb_int<1>> make_ga_vector_fixed(const jtk::vec3<double>& pt, const quantization& q)
    {
    ga::vector<fixint::limb_int<1>> v;
    v.e0 = q.to_int(pt[0]);
    v.e1 = q.to_int(pt[1]);
    v.e2 = q.to_int(pt[2]);
    v.e3 = int(1);
    return v;
    }

  template <int N>
  jtk::vec3<double> to_double3(const ga::vector<typename fixint::bit_int<N>::type>& pt, const quantization& q)
    {
    double x = fixint::to_double(pt.e0);
    double y = fixint::to_double(pt.e1);
    double z = fixint::to_double(pt.e2);
    double w = fixint::to_double(pt.e3);
    x /= w;
    y /= w;
    z /= w;
    jtk::vec3<double> d3;
    d3[0] = q.get_reshrink_factor() * x;
    d3[1] = q.get_reshrink_factor() * y;
    d3[2] = q.get_reshrink_factor() * z;
    return d3;
    }

  bool intersects(bool& degenerate, const jtk::vec2<uint32_t>& e, jtk::vec3<uint32_t> t, const std::vector<jtk::vec3<double>>& vertices, const quantization& q)
    {
    degenerate = false;

    if (t[0] > t[1]) std::swap(t[0], t[1]);
    if (t[1] > t[2]) std::swap(t[1], t[2]);
    if (t[0] > t[1]) std::swap(t[0], t[1]);

    auto ebox = compute_bb(e, vertices);
    auto tbox = compute_bb(t, vertices);

    if (!::intersects(ebox, tbox))
      return false;

    if (has_common_vertex(e, t))
      return false; // not interested in the case where the edge and triangle have a common vertex

    int filter = ga::triangle_edge_intersect_floating_point_filter<double>(make_ga_vector(vertices[t[0]]), make_ga_vector(vertices[t[1]]), make_ga_vector(vertices[t[2]]), make_ga_vector(vertices[e[0]]), make_ga_vector(vertices[e[1]]));
    if (filter == 0) // uncertain, more precision is needed
      {
      filter = triangle_edge_intersect_fixed(make_ga_vector_fixed(vertices[t[0]], q), make_ga_vector_fixed(vertices[t[1]], q), make_ga_vector_fixed(vertices[t[2]], q), make_ga_vector_fixed(vertices[e[0]], q), make_ga_vector_fixed(vertices[e[1]], q));
      if (filter == 0)
        degenerate = true;
      }
    return filter > 0; // intersection
    }

  bool intersects(bool& degenerate, jtk::vec3<uint32_t> t0, jtk::vec3<uint32_t> t1, jtk::vec3<uint32_t> t2, const std::vector<jtk::vec3<double>>& vertices, const quantization& q)
    {
    // This function should only be called if we've already
    // identified that the intersection edges
    //      (t0,t1), (t0,t2), (t1,t2)
    // exist.
    // From this, we can conclude that each pair of triangles
    // shares no more than a single vertex in common.
    //  If each of these shared vertices is different from each other,
    // then we could legitimately have a triple intersection point,
    // but if all three pairs share the same vertex in common, then
    // the intersection of the three triangles must be that vertex.
    // So, we must check for such a single vertex in common amongst
    // the three triangles

    degenerate = false;

    if (t0[0] > t0[1]) std::swap(t0[0], t0[1]);
    if (t0[1] > t0[2]) std::swap(t0[1], t0[2]);
    if (t0[0] > t0[1]) std::swap(t0[0], t0[1]);
    if (t1[0] > t1[1]) std::swap(t1[0], t1[1]);
    if (t1[1] > t1[2]) std::swap(t1[1], t1[2]);
    if (t1[0] > t1[1]) std::swap(t1[0], t1[1]);
    if (t2[0] > t2[1]) std::swap(t2[0], t2[1]);
    if (t2[1] > t2[2]) std::swap(t2[1], t2[2]);
    if (t2[0] > t2[1]) std::swap(t2[0], t2[1]);

    int64_t vcommon = get_common_vertex(t0, t1);
    if (vcommon >= 0)
      {
      for (uint32_t i = 0; i < 3; ++i)
        {
        if (vcommon == t2[i])
          return false; // all three triangles share a common vertex. We're not interested in this case.
        }
      }
    int filter = ga::triangle_triangle_triangle_intersect_floating_point_filter<double>(make_ga_vector(vertices[t0[0]]), make_ga_vector(vertices[t0[1]]), make_ga_vector(vertices[t0[2]]), make_ga_vector(vertices[t1[0]]), make_ga_vector(vertices[t1[1]]), make_ga_vector(vertices[t1[2]]), make_ga_vector(vertices[t2[0]]), make_ga_vector(vertices[t2[1]]), make_ga_vector(vertices[t2[2]]));
    if (filter == 0) // uncertain, more precision is needed
      {
      filter = triangle_triangle_triangle_intersect_fixed(make_ga_vector_fixed(vertices[t0[0]], q), make_ga_vector_fixed(vertices[t0[1]], q), make_ga_vector_fixed(vertices[t0[2]], q), make_ga_vector_fixed(vertices[t1[0]], q), make_ga_vector_fixed(vertices[t1[1]], q), make_ga_vector_fixed(vertices[t1[2]], q), make_ga_vector_fixed(vertices[t2[0]], q), make_ga_vector_fixed(vertices[t2[1]], q), make_ga_vector_fixed(vertices[t2[2]], q));
      if (filter == 0)
        degenerate = true;
      }
    return filter > 0; // intersection
    }

  jtk::vec3<double> compute_intersection(const jtk::vec2<uint32_t>& e, jtk::vec3<uint32_t> t, const std::vector<jtk::vec3<double>>& vertices, const quantization& q)
    {
    if (t[0] > t[1]) std::swap(t[0], t[1]);
    if (t[1] > t[2]) std::swap(t[1], t[2]);
    if (t[0] > t[1]) std::swap(t[0], t[1]);
    auto pt = triangle_edge_intersection_coordinate_fixed(make_ga_vector_fixed(vertices[t[0]], q), make_ga_vector_fixed(vertices[t[1]], q), make_ga_vector_fixed(vertices[t[2]], q), make_ga_vector_fixed(vertices[e[0]], q), make_ga_vector_fixed(vertices[e[1]], q));
    return to_double3<161>(pt, q);
    }

  jtk::vec3<double> compute_intersection(jtk::vec3<uint32_t> t0, jtk::vec3<uint32_t> t1, jtk::vec3<uint32_t> t2, const std::vector<jtk::vec3<double>>& vertices, const quantization& q)
    {
    if (t0[0] > t0[1]) std::swap(t0[0], t0[1]);
    if (t0[1] > t0[2]) std::swap(t0[1], t0[2]);
    if (t0[0] > t0[1]) std::swap(t0[0], t0[1]);
    if (t1[0] > t1[1]) std::swap(t1[0], t1[1]);
    if (t1[1] > t1[2]) std::swap(t1[1], t1[2]);
    if (t1[0] > t1[1]) std::swap(t1[0], t1[1]);
    if (t2[0] > t2[1]) std::swap(t2[0], t2[1]);
    if (t2[1] > t2[2]) std::swap(t2[1], t2[2]);
    if (t2[0] > t2[1]) std::swap(t2[0], t2[1]);
    auto pt = triangle_triangle_triangle_intersection_coordinate_fixed(make_ga_vector_fixed(vertices[t0[0]], q), make_ga_vector_fixed(vertices[t0[1]], q), make_ga_vector_fixed(vertices[t0[2]], q), make_ga_vector_fixed(vertices[t1[0]], q), make_ga_vector_fixed(vertices[t1[1]], q), make_ga_vector_fixed(vertices[t1[2]], q), make_ga_vector_fixed(vertices[t2[0]], q), make_ga_vector_fixed(vertices[t2[1]], q), make_ga_vector_fixed(vertices[t2[2]], q));
    return to_double3<291>(pt, q);
    }

  struct edge_tria_intersecting_pair
    {
    jtk::vec2<uint32_t> e;
    uint32_t t;
    uint32_t v;
    };

  struct triangle_to_intersection
    {
    std::vector<uint32_t> link_to_new_interior_vertices;
    std::vector<std::pair<uint32_t, jtk::vec2<uint32_t>>> link_to_new_boundary_vertices;
    std::vector<edge_tria_intersecting_pair> link_to_edge_tria_intersecting_pair;
    };

  struct triple_intersecting_triangles
    {
    uint32_t t0, t1, t2;
    triple_intersecting_triangles(uint32_t tp0, uint32_t tp1, uint32_t tp2) : t0(tp0), t1(tp1), t2(tp2)
      {}
    };

  std::vector<uint32_t> get_unique_triangles(const std::vector<edge_tria_intersecting_pair>& intersections)
    {
    std::vector<uint32_t> tria;
    tria.reserve(intersections.size());
    for (const auto& in : intersections)
      tria.push_back(in.t);
    std::sort(tria.begin(), tria.end());
    tria.erase(std::unique(tria.begin(), tria.end()), tria.end());
    return tria;
    }

  bool _compute_intersection_data(std::vector<std::vector<std::pair<jtk::vec2<uint32_t>, uint32_t>>>& intersections_per_tria, const std::vector<triangle_info>& info, const e_t_map& edge_map, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<double>>& vertices, bool resolve_all_intersections, const cork_options& options, const quantization& q)
    {
    FILE* f = nullptr;
    if (options.debug_folder)
      {
      std::string filename = _get_debug_file("find_intersections.txt", options);
      if (options.p_str)
        *options.p_str << "[LIBCORK] " << "Writing _compute_intersection_data debug output to " << filename << "\n";
      f = fopen(filename.c_str(), "w");
      }
    intersections_per_tria.clear();
    intersections_per_tria.resize(triangles.size());
    std::atomic<bool> degen{ false };
    jtk::timer t;
    t.start();
    int runs = resolve_all_intersections ? 1 : 2;
    for (int m = 0; m < runs; ++m)
      {
      t.start();
      int mesh_id = resolve_all_intersections ? -1 : m;
      auto edge_aabvh = compute_edge_aabvh(edge_map, vertices, info, mesh_id);
      double time_elapsed = t.time_elapsed();
      if (options.p_str)
        *options.p_str << "[LIBCORK] " << "Computing edge aabvh took " << time_elapsed << "s\n";

      if (f)
        {
        fprintf(f, "Creating edge BVH took %fs\n", time_elapsed);
        fprintf(f, "The resulting edge BVH is has %d leaves\n", (int)edge_aabvh.number_of_leafs());
        }

      t.start();

      aabvh_edge_triangle_parallel([&](const jtk::vec2<uint32_t>& e, const jtk::vec3<uint32_t>& t, uint32_t tria_index)->bool
        {
        bool degenerate;
        if (intersects(degenerate, e, t, vertices, q))
          {
          intersections_per_tria[tria_index].push_back(std::pair<jtk::vec2<uint32_t>, uint32_t>(e, tria_index));
          }
        if (degenerate)
          {
          degen = true;
          return false;
          }
        return true; // continue with aabvh_edge_triangle
        }, edge_aabvh, info, triangles, vertices, mesh_id, options);

      time_elapsed = t.time_elapsed();

      if (options.p_str)
        *options.p_str << "[LIBCORK] " << "Finding all edge-triangle intersections took " << time_elapsed << "s\n";
      if (f)
        fprintf(f, "Finding all edge-triangle intersections took %fs\n", time_elapsed);
      if (degen)
        break;
      }
    if (f)
      {
      uint32_t intersections_with_mesh_0 = 0;
      uint32_t intersections_with_mesh_1 = 0;
      for (uint32_t tr = 0; tr < (uint32_t)triangles.size(); ++tr)
        {
        const auto& isct = intersections_per_tria[tr];
        fprintf(f, "Triangle %d (mesh %d) has %d intersections: ", (int)tr, info[tr].bit & 1, (int)isct.size());
        if ((info[tr].bit & 1) == 0)
          intersections_with_mesh_0 += (uint32_t)isct.size();
        else
          intersections_with_mesh_1 += (uint32_t)isct.size();
        for (auto e : isct)
          fprintf(f, "(%d, %d) ", (int)e.first[0], (int)e.first[1]);
        fprintf(f, "\n");
        }
      fprintf(f, "Found %d intersections with first mesh\n", (int)intersections_with_mesh_0);
      fprintf(f, "Found %d intersections with second mesh\n", (int)intersections_with_mesh_1);
      }

    if (f)
      fclose(f);
    return !degen;
    }

  bool try_to_find_intersections(std::vector<triangle_to_intersection>& intersection_data, std::vector<jtk::vec3<double>>& new_vertices, const std::vector<triangle_info>& info, const e_t_map& edge_map, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<double>>& vertices, bool resolve_all_intersections, const cork_options& options, const quantization& q)
    {
    intersection_data.clear();
    new_vertices.clear();
    intersection_data.resize(triangles.size());

    std::vector<std::vector<std::pair<jtk::vec2<uint32_t>, uint32_t>>> intersections_per_tria;
    if (!_compute_intersection_data(intersections_per_tria, info, edge_map, triangles, vertices, resolve_all_intersections, options, q))
      {
      return false;
      }

    std::vector<std::pair<jtk::vec2<uint32_t>, uint32_t>> intersections;
    for (const auto& intersect : intersections_per_tria)
      intersections.insert(intersections.end(), intersect.begin(), intersect.end());
    if (options.p_str)
      *options.p_str << "[LIBCORK] " << "Found " << (int)intersections.size() << " intersections\n";

    new_vertices.resize(intersections.size());
    jtk::timer t;
    t.start();
    _parallel_for((uint32_t)0, (uint32_t)intersections.size(), [&](uint32_t new_vertex_id)
      {
      auto pr = intersections[new_vertex_id];
      auto e = pr.first;
      uint32_t tria_index = pr.second;
      auto t = triangles[tria_index];
      new_vertices[new_vertex_id] = compute_intersection(e, t, vertices, q);
      }, options);
    for (uint32_t new_vertex_id = 0; new_vertex_id < (uint32_t)intersections.size(); ++new_vertex_id)
      {
      auto pr = intersections[new_vertex_id];
      auto e = pr.first;
      uint32_t tria_index = pr.second;
      intersection_data[tria_index].link_to_new_interior_vertices.push_back(new_vertex_id);
      std::vector<uint32_t> edge_trias = edge_map[e[0]].find(e[1])->second;
      for (auto tri : edge_trias)
        {
        edge_tria_intersecting_pair etp;
        etp.e = e;
        etp.t = tri;
        etp.v = new_vertex_id;
        intersection_data[tria_index].link_to_edge_tria_intersecting_pair.push_back(etp);
        intersection_data[tri].link_to_new_boundary_vertices.push_back(std::pair<uint32_t, jtk::vec2<uint32_t>>(new_vertex_id, e));
        etp.t = tria_index;
        intersection_data[tri].link_to_edge_tria_intersecting_pair.push_back(etp);
        }
      }
    if (options.p_str)
      *options.p_str << "[LIBCORK] " << "Computing all edge-triangle intersections took " << t.time_elapsed() << "s\n";

    t.start();
    std::vector<triple_intersecting_triangles> triples;
    for (uint32_t t0 = 0; t0 < (uint32_t)intersection_data.size(); ++t0)
      {
      const auto& isct0 = intersection_data[t0];
      auto isct0_triangles = get_unique_triangles(isct0.link_to_edge_tria_intersecting_pair);
      for (uint32_t i = 0; i < (uint32_t)isct0_triangles.size(); ++i)
        for (uint32_t j = i + 1; j < (uint32_t)isct0_triangles.size(); ++j)
          {
          uint32_t t1 = isct0_triangles[i];
          uint32_t t2 = isct0_triangles[j];
          if (t0 < t1 && t0 < t2)
            {
            const auto& isct1 = intersection_data[t1];
            auto isct1_triangles = get_unique_triangles(isct1.link_to_edge_tria_intersecting_pair);
            for (auto tr : isct1_triangles)
              {
              if (tr == t2)
                {
                triples.push_back(triple_intersecting_triangles(t0, t1, t2));
                }
              }
            }
          }
      }
    if (options.p_str)
      *options.p_str << "[LIBCORK] " << "Finding all triangle-triangle-triangle intersections took " << t.time_elapsed() << "s\n";

    t.start();
    int triple_intersections = 0;
    for (auto triple : triples)
      {
      bool degen;
      if (!intersects(degen, triangles[triple.t0], triangles[triple.t1], triangles[triple.t2], vertices, q))
        {
        if (degen)
          return false;
        continue;
        }
      ++triple_intersections;
      uint32_t new_vertex_id = (uint32_t)new_vertices.size();
      new_vertices.push_back(compute_intersection(triangles[triple.t0], triangles[triple.t1], triangles[triple.t2], vertices, q));
      intersection_data[triple.t0].link_to_new_interior_vertices.push_back(new_vertex_id);
      intersection_data[triple.t1].link_to_new_interior_vertices.push_back(new_vertex_id);
      intersection_data[triple.t2].link_to_new_interior_vertices.push_back(new_vertex_id);
      edge_tria_intersecting_pair etp;
      etp.t = triple.t1;
      etp.v = new_vertex_id;
      intersection_data[triple.t0].link_to_edge_tria_intersecting_pair.push_back(etp);
      intersection_data[triple.t2].link_to_edge_tria_intersecting_pair.push_back(etp);
      etp.t = triple.t2;
      intersection_data[triple.t0].link_to_edge_tria_intersecting_pair.push_back(etp);
      intersection_data[triple.t1].link_to_edge_tria_intersecting_pair.push_back(etp);
      etp.t = triple.t0;
      intersection_data[triple.t2].link_to_edge_tria_intersecting_pair.push_back(etp);
      intersection_data[triple.t1].link_to_edge_tria_intersecting_pair.push_back(etp);
      }
    if (options.p_str)
      {
      *options.p_str << "[LIBCORK] " << "Found " << triple_intersections << " triple intersections\n";
      *options.p_str << "[LIBCORK] " << "Computing all triangle-triangle-triangle intersections took " << t.time_elapsed() << "s\n";
      }

    return true;
    }

  void _subdivide_triangles(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<triangle_info>& info, std::vector<jtk::vec3<double>>& vertices, const std::vector<triangle_to_intersection>& intersection_data, std::vector<jtk::vec3<double>>& new_vertices, const cork_options& options)
    {
    FILE* f = nullptr;
    if (options.debug_folder)
      {
      std::string filename = _get_debug_file("subdivide.txt", options);
      if (options.p_str)
        *options.p_str << "[LIBCORK] " << "Writing _subdivide_triangles debug output to " << filename << "\n";
      f = fopen(filename.c_str(), "w");
      }

    std::mutex mut;

    uint32_t original_nr_of_vertices = (uint32_t)vertices.size();
    vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
    uint32_t original_nr_of_triangles = (uint32_t)triangles.size();
    std::vector<std::vector<jtk::vec3<uint32_t>>> new_triangles(triangles.size());
    std::vector<jtk::vec3<double>> vertices_from_segment_intersections;
    uint32_t new_vertex_id = (uint32_t)vertices.size();
    _parallel_for((uint32_t)0, original_nr_of_triangles, [&](uint32_t t)
      {
      const jtk::vec3<uint32_t>& tria = triangles[t];
      const triangle_to_intersection& isct = intersection_data[t];
      if (isct.link_to_new_boundary_vertices.empty() && isct.link_to_new_interior_vertices.empty())
        return;
      if (f)
        {
        fprintf(f, "Working on triangle %d\n", (int)t);
        fprintf(f, "V0: (%f, %f, %f)\n", vertices[tria[0]][0], vertices[tria[0]][1], vertices[tria[0]][2]);
        fprintf(f, "V1: (%f, %f, %f)\n", vertices[tria[1]][0], vertices[tria[1]][1], vertices[tria[1]][2]);
        fprintf(f, "V2: (%f, %f, %f)\n", vertices[tria[2]][0], vertices[tria[2]][1], vertices[tria[2]][2]);
        }
      for (auto iter = isct.link_to_new_interior_vertices.begin(); iter != isct.link_to_new_interior_vertices.end(); ++iter)
        {
        auto l = *iter;
        if (f)
          fprintf(f, "Interior vertex: (%f, %f, %f) interior\n", new_vertices[l][0], new_vertices[l][1], new_vertices[l][2]);
        }
      for (auto iter = isct.link_to_new_boundary_vertices.rbegin(); iter != isct.link_to_new_boundary_vertices.rend(); ++iter)
        {
        auto l = *iter;
        if (f)
          fprintf(f, "Interior vertex: (%f, %f, %f) boundary\n", new_vertices[l.first][0], new_vertices[l.first][1], new_vertices[l.first][2]);
        }
      jtk::vec2<uint32_t> e0(tria[0], tria[1]);
      if (e0[1] < e0[0])
        std::swap(e0[1], e0[0]);
      jtk::vec2<uint32_t> e1(tria[1], tria[2]);
      if (e1[1] < e1[0])
        std::swap(e1[1], e1[0]);
      jtk::vec2<uint32_t> e2(tria[2], tria[0]);
      if (e2[1] < e2[0])
        std::swap(e2[1], e2[0]);
      std::vector<int> boundary_vertices_on_edge_0;
      std::vector<int> boundary_vertices_on_edge_1;
      std::vector<int> boundary_vertices_on_edge_2;
      std::vector<int> local_vertices;
      local_vertices.push_back(tria[0]);
      local_vertices.push_back(tria[1]);
      local_vertices.push_back(tria[2]);
      boundary_vertices_on_edge_0.push_back(0);
      boundary_vertices_on_edge_0.push_back(1);
      boundary_vertices_on_edge_1.push_back(1);
      boundary_vertices_on_edge_1.push_back(2);
      boundary_vertices_on_edge_2.push_back(2);
      boundary_vertices_on_edge_2.push_back(0);
      int index = 3;
      for (auto pr : isct.link_to_new_boundary_vertices)
        {
        if (pr.second == e0)
          {
          local_vertices.push_back(pr.first + original_nr_of_vertices);
          boundary_vertices_on_edge_0.push_back(index++);
          }
        else if (pr.second == e1)
          {
          local_vertices.push_back(pr.first + original_nr_of_vertices);
          boundary_vertices_on_edge_1.push_back(index++);
          }
        else if (pr.second == e2)
          {
          local_vertices.push_back(pr.first + original_nr_of_vertices);
          boundary_vertices_on_edge_2.push_back(index++);
          }
        else
          {
          if (options.p_str)
            *options.p_str << "Error: could not locate boundary vertex" << std::endl;
          exit(1);
          }
        }
      for (auto v : isct.link_to_new_interior_vertices)
        local_vertices.push_back(v + original_nr_of_vertices);

      std::vector<jtk::vec2<double>> projected_points;

      auto normal = jtk::cross(vertices[tria[1]] - vertices[tria[0]], vertices[tria[2]] - vertices[tria[0]]);
      int norm_dim = 0;
      if (std::abs(normal[1]) > std::abs(normal[norm_dim]))
        norm_dim = 1;
      if (std::abs(normal[2]) > std::abs(normal[norm_dim]))
        norm_dim = 2;
      int dim0 = (norm_dim + 1) % 3;
      int dim1 = (norm_dim + 2) % 3;
      double sign_flip = (normal[norm_dim] < 0.0) ? -1.0 : 1.0;

      std::map<uint32_t, uint32_t> global_to_local;
      for (uint32_t ind = 0; ind < (uint32_t)local_vertices.size(); ++ind)
        {
        auto v = local_vertices[ind];
        projected_points.push_back(jtk::vec2<double>(vertices[v][dim0], vertices[v][dim1] * sign_flip));
        global_to_local[v] = ind;
        }
      std::vector<jtk::vec2<uint32_t>> interior_edges;

      std::map<uint32_t, uint32_t> tria_to_index;
      std::vector<std::vector<uint32_t>> connected_vertices;
      for (uint32_t i = 0; i < (uint32_t)isct.link_to_edge_tria_intersecting_pair.size(); ++i)
        {
        for (uint32_t j = i + 1; j < (uint32_t)isct.link_to_edge_tria_intersecting_pair.size(); ++j)
          {
          uint32_t t1 = isct.link_to_edge_tria_intersecting_pair[i].t;
          uint32_t t2 = isct.link_to_edge_tria_intersecting_pair[j].t;
          if (t1 == t2) // create interior edge
            {
            auto it = tria_to_index.find(t1);
            if (it == tria_to_index.end()) // first encounter of this triangle
              {
              uint32_t vindex = (uint32_t)connected_vertices.size();
              tria_to_index[t1] = vindex;
              uint32_t v1 = isct.link_to_edge_tria_intersecting_pair[i].v + original_nr_of_vertices;
              uint32_t v2 = isct.link_to_edge_tria_intersecting_pair[j].v + original_nr_of_vertices;
              connected_vertices.emplace_back();
              connected_vertices.back().push_back(v1);
              connected_vertices.back().push_back(v2);
              }
            else
              {
              uint32_t v1 = isct.link_to_edge_tria_intersecting_pair[i].v + original_nr_of_vertices;
              uint32_t v2 = isct.link_to_edge_tria_intersecting_pair[j].v + original_nr_of_vertices;
              connected_vertices[it->second].push_back(v1);
              connected_vertices[it->second].push_back(v2);
              }
            }
          }
        }

      for (auto v_list : connected_vertices)
        {
        if (v_list.size() == 2)
          {
          uint32_t v1 = v_list[0];
          uint32_t v2 = v_list[1];
          if (global_to_local.find(v1) == global_to_local.end())
            {
            if (options.p_str)
              *options.p_str << "[LIBCORK] ERROR" << std::endl;
            }
          if (global_to_local.find(v2) == global_to_local.end())
            {
            if (options.p_str)
              *options.p_str << "[LIBCORK] ERROR" << std::endl;
            }
          jtk::vec2<uint32_t> e(global_to_local[v1], global_to_local[v2]);
          interior_edges.push_back(e);
          }
        else
          {
          std::sort(v_list.begin(), v_list.end());
          v_list.erase(std::unique(v_list.begin(), v_list.end()), v_list.end());
          std::vector<uint32_t> local_ids;
          for (auto v : v_list)
            {
            if (global_to_local.find(v) == global_to_local.end())
              {
              if (options.p_str)
                *options.p_str << "[LIBCORK] ERROR" << std::endl;
              }
            local_ids.push_back(global_to_local[v]);
            }
          auto V0 = projected_points[local_ids[0]];
          auto V1 = projected_points[local_ids[1]];
          auto V1V0 = V1 - V0;
          int alpha_dim = 0;
          if (std::abs(V1V0[1]) > std::abs(V1V0[0]))
            alpha_dim = 1;
          std::vector<double> alpha;
          alpha.push_back(0.0);
          alpha.push_back(1.0);
          for (uint32_t i = 2; i < (uint32_t)local_ids.size(); ++i)
            {
            auto P = projected_points[local_ids[i]];
            double cur_alpha = (P[alpha_dim] - V0[alpha_dim]) / V1V0[alpha_dim];
            alpha.push_back(cur_alpha);
            }
          auto segment_order = sort_indices(alpha);
          for (uint32_t i = 0; i < (uint32_t)segment_order.size() - 1; ++i)
            {
            jtk::vec2<uint32_t> e(local_ids[segment_order[i]], local_ids[segment_order[i + 1]]);
            interior_edges.push_back(e);
            }
          }
        }

      struct triangulateio in, out;

      in.numberofpoints = (int)projected_points.size();
      in.numberofpointattributes = 0;
      in.pointlist = new REAL[in.numberofpoints * 2];
      in.pointattributelist = nullptr;
      in.pointmarkerlist = new int[in.numberofpoints];
      for (int i = 0; i < in.numberofpoints; ++i)
        {
        in.pointlist[i * 2 + 0] = projected_points[i][0];
        in.pointlist[i * 2 + 1] = projected_points[i][1];
        in.pointmarkerlist[i] = 0;
        }
      for (auto v : boundary_vertices_on_edge_0)
        in.pointmarkerlist[v] = 1;
      for (auto v : boundary_vertices_on_edge_1)
        in.pointmarkerlist[v] = 1;
      for (auto v : boundary_vertices_on_edge_2)
        in.pointmarkerlist[v] = 1;
      if (f)
        {
        for (int k = 0; k < in.numberofpoints; ++k)
          {
          fprintf(f, "Delaunay pt (%f, %f) - b:%d\n", in.pointlist[k * 2], in.pointlist[k * 2 + 1], in.pointmarkerlist[k]);
          }
        }
      in.numberofsegments = (int)(boundary_vertices_on_edge_0.size() + boundary_vertices_on_edge_1.size() + boundary_vertices_on_edge_2.size() - 3 + interior_edges.size());
      in.numberofholes = 0;
      in.numberofregions = 0;
      in.segmentlist = new int[in.numberofsegments * 2];
      in.segmentmarkerlist = new int[in.numberofsegments];
      for (int i = 0; i < in.numberofsegments; ++i)
        in.segmentmarkerlist[i] = 1;

      int alpha_dim;
      jtk::vec2<double> V0, V1V0;

      std::vector<double> boundary_vertices_on_edge_0_alpha(boundary_vertices_on_edge_0.size(), 0.0);
      boundary_vertices_on_edge_0_alpha[1] = 1.0;
      V0 = projected_points[boundary_vertices_on_edge_0[0]];
      V1V0 = projected_points[boundary_vertices_on_edge_0[1]] - V0;
      alpha_dim = 0;
      if (std::abs(V1V0[1]) > std::abs(V1V0[0]))
        alpha_dim = 1;
      for (int i = 2; i < boundary_vertices_on_edge_0.size(); ++i)
        {
        auto P = projected_points[boundary_vertices_on_edge_0[i]];
        double alpha = (P[alpha_dim] - V0[alpha_dim]) / V1V0[alpha_dim];
        if (alpha <= 0.0 || alpha >= 1.0 || std::find(boundary_vertices_on_edge_0_alpha.begin(), boundary_vertices_on_edge_0_alpha.end(), alpha) != boundary_vertices_on_edge_0_alpha.end())
          {
          if (options.p_str)
            *options.p_str << "[LIBCORK] ERROR: alpha value already exists" << std::endl;
          if (f)
            fprintf(f, "ERROR: alpha value already exists\n");
          }
        boundary_vertices_on_edge_0_alpha[i] = alpha;
        }
      std::vector<double> boundary_vertices_on_edge_1_alpha(boundary_vertices_on_edge_1.size(), 0.0);
      boundary_vertices_on_edge_1_alpha[1] = 1.0;
      V0 = projected_points[boundary_vertices_on_edge_1[0]];
      V1V0 = projected_points[boundary_vertices_on_edge_1[1]] - V0;
      alpha_dim = 0;
      if (std::abs(V1V0[1]) > std::abs(V1V0[0]))
        alpha_dim = 1;
      for (int i = 2; i < boundary_vertices_on_edge_1.size(); ++i)
        {
        auto P = projected_points[boundary_vertices_on_edge_1[i]];
        double alpha = (P[alpha_dim] - V0[alpha_dim]) / V1V0[alpha_dim];
        if (alpha <= 0.0 || alpha >= 1.0 || std::find(boundary_vertices_on_edge_1.begin(), boundary_vertices_on_edge_1.end(), alpha) != boundary_vertices_on_edge_1.end())
          {
          if (options.p_str)
            *options.p_str << "[LIBCORK] ERROR: alpha value already exists" << std::endl;
          if (f)
            fprintf(f, "ERROR: alpha value already exists\n");
          }
        boundary_vertices_on_edge_1_alpha[i] = alpha;
        }
      std::vector<double> boundary_vertices_on_edge_2_alpha(boundary_vertices_on_edge_2.size(), 0.0);
      boundary_vertices_on_edge_2_alpha[1] = 1.0;
      V0 = projected_points[boundary_vertices_on_edge_2[0]];
      V1V0 = projected_points[boundary_vertices_on_edge_2[1]] - V0;
      alpha_dim = 0;
      if (std::abs(V1V0[1]) > std::abs(V1V0[0]))
        alpha_dim = 1;
      for (int i = 2; i < boundary_vertices_on_edge_2.size(); ++i)
        {
        auto P = projected_points[boundary_vertices_on_edge_2[i]];
        double alpha = (P[alpha_dim] - V0[alpha_dim]) / V1V0[alpha_dim];
        if (alpha <= 0.0 || alpha >= 1.0 || std::find(boundary_vertices_on_edge_2_alpha.begin(), boundary_vertices_on_edge_2_alpha.end(), alpha) != boundary_vertices_on_edge_2_alpha.end())
          {
          if (options.p_str)
            *options.p_str << "[LIBCORK] ERROR: alpha value already exists" << std::endl;
          if (f)
            fprintf(f, "ERROR: alpha value already exists\n");
          }
        boundary_vertices_on_edge_2_alpha[i] = alpha;
        }

      auto segment_order_edge0 = sort_indices(boundary_vertices_on_edge_0_alpha);
      auto segment_order_edge1 = sort_indices(boundary_vertices_on_edge_1_alpha);
      auto segment_order_edge2 = sort_indices(boundary_vertices_on_edge_2_alpha);

      index = 0;
      for (int i = 0; i < segment_order_edge0.size() - 1; ++i)
        {
        uint32_t v0 = boundary_vertices_on_edge_0[segment_order_edge0[i]];
        uint32_t v1 = boundary_vertices_on_edge_0[segment_order_edge0[i + 1]];
        in.segmentlist[index * 2] = v0;
        in.segmentlist[index * 2 + 1] = v1;
        ++index;
        }
      for (int i = 0; i < segment_order_edge1.size() - 1; ++i)
        {
        uint32_t v0 = boundary_vertices_on_edge_1[segment_order_edge1[i]];
        uint32_t v1 = boundary_vertices_on_edge_1[segment_order_edge1[i + 1]];
        in.segmentlist[index * 2] = v0;
        in.segmentlist[index * 2 + 1] = v1;
        ++index;
        }
      for (int i = 0; i < segment_order_edge2.size() - 1; ++i)
        {
        uint32_t v0 = boundary_vertices_on_edge_2[segment_order_edge2[i]];
        uint32_t v1 = boundary_vertices_on_edge_2[segment_order_edge2[i + 1]];
        in.segmentlist[index * 2] = v0;
        in.segmentlist[index * 2 + 1] = v1;
        ++index;
        }
      for (int i = 0; i < interior_edges.size(); ++i)
        {
        in.segmentlist[index * 2] = interior_edges[i][0];
        in.segmentlist[index * 2 + 1] = interior_edges[i][1];
        in.segmentmarkerlist[index] = 0;
        ++index;
        }

      if (f)
        {
        for (int i = 0; i < in.numberofsegments; ++i)
          {
          fprintf(f, "Segment (%d, %d) b:%d\n", in.segmentlist[i * 2], in.segmentlist[i * 2 + 1], in.segmentmarkerlist[i]);
          }
        }
      in.numberoftriangles = 0;
      in.numberoftriangleattributes = 0;

      out.pointlist = nullptr;
      out.pointattributelist = nullptr;
      out.pointmarkerlist = nullptr;
      out.trianglelist = nullptr;
      out.segmentlist = nullptr;
      out.segmentmarkerlist = nullptr;


      char* params = (char*)("pzQYY");
      //char* debug_params = (char*)("pzYYVC");
      triangulate(params, &in, &out, nullptr);

      if (out.numberofpoints < in.numberofpoints)
        {
        if (options.p_str)
          {
          *options.p_str << "[LIBCORK] ERROR: Could not triangulate all the input points\n";
          }
        if (f)
          fprintf(f, "ERROR: Could not triangulate all the input points\n");
        exit(1);
        }

      if (out.numberofpoints > in.numberofpoints)
        {
        if (options.p_str)
          {
          *options.p_str << "[LIBCORK] Subdivision created extra vertices for triangle " << t << std::endl;
          }
        if (f)
          {
          fprintf(f, "out.numberofpoints: %d\n", out.numberofpoints);
          fprintf(f, "projected_points.size(): %d\n", (int)projected_points.size());
          fprintf(f, "dumping out the points' coordinates\n");
          for (int k = 0; k < projected_points.size(); ++k)
            fprintf(f, "%d: %f %f\n", k, projected_points[k][0], projected_points[k][1]);
          fprintf(f, "dumping out the segments\n");
          for (int k = 0; k < in.numberofsegments; ++k)
            fprintf(f, "  %d; %d (%d)\n", in.segmentlist[k * 2 + 0], in.segmentlist[k * 2 + 1], in.segmentmarkerlist[k]);
          fprintf(f, "dumping out the solved for triangles\n");
          for (int k = 0; k < out.numberoftriangles; ++k)
            fprintf(f, " %d; %d; %d\n", out.trianglelist[(k * 3) + 0], out.trianglelist[(k * 3) + 1], out.trianglelist[(k * 3) + 2]);
          }
        for (int i = 0; i < in.numberofpoints; ++i)
          {
          if ((projected_points[i][0] != out.pointlist[i * 2]) || (projected_points[i][1] != out.pointlist[i * 2 + 1]))
            {
            printf("ERROR: Point moved\n");
            if (f)
              fprintf(f, "ERROR: Point moved\n");
            }
          }
        for (int i = in.numberofpoints; i < out.numberofpoints; ++i)
          {
          // Compute barycentric coordinates (u,v,w) for point p with respect to triangle (a,b,c)
          jtk::vec2<double> a = projected_points[0];
          jtk::vec2<double> b = projected_points[1];
          jtk::vec2<double> c = projected_points[2];
          jtk::vec2<double> p = jtk::vec2<double>(out.pointlist[i * 2], out.pointlist[i * 2 + 1]);
          auto v0 = b - a;
          auto v1 = c - a;
          auto v2 = p - a;
          double den = v0[0] * v1[1] - v1[0] * v0[1];
          double v = (v2[0] * v1[1] - v1[0] * v2[1]) / den;
          double w = (v0[0] * v2[1] - v2[0] * v0[1]) / den;
          double u = 1.0 - v - w;
          if (options.p_str)
            *options.p_str << "[LIBCORK] Barycentric coordinates are (" << u << ", " << v << ", " << w << ")" << std::endl;
          jtk::vec3<double> pt = u * vertices[tria[0]] + v * vertices[tria[1]] + w * vertices[tria[2]];
          mut.lock();
          local_vertices.push_back(new_vertex_id);
          ++new_vertex_id;
          vertices_from_segment_intersections.push_back(pt);
          info[t].marked = true;
          for (auto pr : isct.link_to_edge_tria_intersecting_pair)
            info[pr.t].marked = true;
          mut.unlock();
          }
        }

      for (int i = 0; i < in.numberofpoints; ++i)
        {
        if (in.pointmarkerlist[i] != out.pointmarkerlist[i])
          {
          if (options.p_str)
            *options.p_str << "[LIBCORK] Created a boundary vertex" << std::endl;
          fprintf(f, "ERROR: created a boundary vertex\n");
          }
        }

      for (int k = 0; k < out.numberoftriangles; ++k)
        {
        int v0 = out.trianglelist[(k * 3) + 0];
        int v1 = out.trianglelist[(k * 3) + 1];
        int v2 = out.trianglelist[(k * 3) + 2];
        jtk::vec3<uint32_t> new_t = jtk::vec3<uint32_t>(local_vertices[v0], local_vertices[v1], local_vertices[v2]);
        new_triangles[t].push_back(new_t);
        }

      free(in.pointlist);
      free(in.pointmarkerlist);
      free(in.segmentlist);
      free(in.segmentmarkerlist);
      free(out.pointlist);
      free(out.pointmarkerlist);
      free(out.trianglelist);
      free(out.segmentlist);
      free(out.segmentmarkerlist);     

      }, options);

    vertices.insert(vertices.end(), vertices_from_segment_intersections.begin(), vertices_from_segment_intersections.end());
    new_vertices.insert(new_vertices.end(), vertices_from_segment_intersections.begin(), vertices_from_segment_intersections.end());

    for (uint32_t t = 0; t < (uint32_t)new_triangles.size(); ++t)
      {
      const auto& new_tri = new_triangles[t];
      if (new_tri.empty())
        continue;
      auto& tria = triangles[t];
      tria = new_tri.front();
      for (int k = 1; k < new_tri.size(); ++k)
        {
        triangles.push_back(new_tri[k]);
        info.push_back(info[t]);
        }
      }

    if (f)
      fclose(f);

    }

  void _resolve_intersections(std::vector<triangle_info>& info, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<double>>& vertices, bool resolve_all_intersections, const cork_options& options)
    {
    jtk::vec3<double> center(0, 0, 0);
    if (options.center)
      {
      for (const auto& pt : vertices)
        center = center + pt;
      center = center / (double)vertices.size();
      for (auto& pt : vertices)
        {
        pt = pt - center;
        }
      }

    srand(100);
    if (info.empty())
      info.resize(triangles.size());
    quantization q(get_magnitude(vertices), 30);
    auto quantized = get_scaled_vertices(vertices, q);
    auto edge_to_tria_map = edges_to_triangles_map((uint32_t)vertices.size(), triangles);

    std::vector<triangle_to_intersection> intersection_data;
    std::vector<jtk::vec3<double>> new_vertices;

    int trys;
    for (trys = 0; trys < 5; ++trys)
      {
      if (try_to_find_intersections(intersection_data, new_vertices, info, edge_to_tria_map, triangles, quantized, resolve_all_intersections, options, q))
        break;
      if (options.p_str)
        *options.p_str << "[LIBCORK] " << "Perturbing vertices to deal with degeneracies" << std::endl;
      perturb_vertices(quantized, q);
      }
    if (trys == 5)
      {
      if (options.p_str)
        *options.p_str << "[LIBCORK] Error: could not solve intersections by perturbing the mesh" << std::endl;
      return;
      }

    jtk::timer t;
    t.start();
    _subdivide_triangles(triangles, info, quantized, intersection_data, new_vertices, options);
    assert(triangles.size() == info.size());
    if (options.p_str)
      *options.p_str << "[LIBCORK] " << "Subdividing triangles took " << t.time_elapsed() << "s\n";
    //vertices.insert(vertices.end(), new_vertices.begin(), new_vertices.end());
    vertices = quantized;
    if (options.center)
      {
      for (auto& pt : vertices)
        {
        pt = pt + center;
        }
      }
    }

  inline jtk::vec3<float> mat_vec_multiply(const float* m, const jtk::vec3<float>& v)
    {
    jtk::vec3<float> out;
    out[0] = m[0] * v[0] + m[4] * v[1] + m[8] * v[2] + m[12];
    out[1] = m[1] * v[0] + m[5] * v[1] + m[9] * v[2] + m[13];
    out[2] = m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14];
    return out;
    }

  void _resolve_intersections(std::vector<jtk::vec3<double>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<triangle_info>& info, const std::vector<jtk::vec3<uint32_t>>& triangles_left, const std::vector<jtk::vec3<float>>& vertices_left, const float* cs_left, const std::vector<jtk::vec3<uint32_t>>& triangles_right, const std::vector<jtk::vec3<float>>& vertices_right, const float* cs_right, const cork_options& options)
    {
    vertices.clear();
    vertices.reserve(vertices_left.size() + vertices_right.size());
    for (const auto& pt : vertices_left)
      {
      auto pttr = mat_vec_multiply(cs_left, pt);
      jtk::vec3<double> v((double)pttr[0], (double)pttr[1], (double)pttr[2]);
      vertices.push_back(v);
      }
    for (const auto& pt : vertices_right)
      {
      auto pttr = mat_vec_multiply(cs_right, pt);
      jtk::vec3<double> v((double)pttr[0], (double)pttr[1], (double)pttr[2]);
      vertices.push_back(v);
      }
    triangles = triangles_left;
    triangles.reserve(triangles_left.size() + triangles_right.size());
    const uint32_t offset = (uint32_t)vertices_left.size();
    for (auto t : triangles_right)
      {
      t[0] += offset;
      t[1] += offset;
      t[2] += offset;
      triangles.push_back(t);
      }

    info.resize(triangles.size());
    for (int t = 0; t < triangles_left.size(); ++t)
      {
      info[t].parent = t;
      info[t].bit = 0;
      info[t].marked = false;
      }
    for (int t = 0; t < triangles_right.size(); ++t)
      {
      info[t + triangles_left.size()].parent = t;
      info[t + triangles_left.size()].bit = 1;
      info[t + triangles_left.size()].marked = false;
      }
    _resolve_intersections(info, triangles, vertices, options.resolve_all_intersections, options);
    }

  bool is_edge_on_intersection(const jtk::vec2<uint32_t>& e, const std::vector<triangle_info>& info, const e_t_map& edge_map)
    {
    auto e0 = e[0];
    auto e1 = e[1];
    if (e1 < e0)
      std::swap(e1, e0);
    assert(edge_map[e0].find(e1) != edge_map[e0].end());
    const auto& trias = edge_map[e0].find(e1)->second;
    assert(!trias.empty());
    auto orig_mesh = info[trias.front()].bit & 1;
    for (uint32_t j = 1; j < (uint32_t)trias.size(); ++j)
      {
      if ((info[trias[j]].bit & 1) != orig_mesh)
        return true;
      }
    return false;
    }

  bool is_edge_with_marked_triangle(const jtk::vec2<uint32_t>& e, const std::vector<triangle_info>& info, const e_t_map& edge_map)
    {
    auto e0 = e[0];
    auto e1 = e[1];
    if (e1 < e0)
      std::swap(e1, e0);
    assert(edge_map[e0].find(e1) != edge_map[e0].end());
    const auto& trias = edge_map[e0].find(e1)->second;
    assert(!trias.empty());
    for (uint32_t j = 0; j < (uint32_t)trias.size(); ++j)
      {
      if (info[trias[j]].marked)
        return true;
      }
    return false;
    }

  bool is_edge_manifold(const jtk::vec2<uint32_t>& e, const std::vector<triangle_info>& info, const e_t_map& edge_map, uint8_t mesh_id)
    {
    auto e0 = e[0];
    auto e1 = e[1];
    if (e1 < e0)
      std::swap(e1, e0);
    assert(edge_map[e0].find(e1) != edge_map[e0].end());
    const auto& trias = edge_map[e0].find(e1)->second;
    assert(!trias.empty());
    int valid_trias = 0;
    for (uint32_t j = 0; j < (uint32_t)trias.size(); ++j)
      {
      if ((info[trias[j]].bit & 1) == mesh_id)
        ++valid_trias;
      }
    return valid_trias == 2;
    }

  std::vector<std::vector<uint32_t>> _compute_components(std::vector<triangle_info>& info, const std::vector<jtk::vec3<uint32_t>>& triangles, const e_t_map& edge_map)
    {
    union_find ds((uint32_t)triangles.size());
    for (uint32_t i = 0; i < (uint32_t)edge_map.size(); ++i)
      {
      for (auto it = edge_map[i].begin(); it != edge_map[i].end(); ++it)
        {
        int tid0 = int(it->second[0]);
        int tid0_ = -1;
        auto bit = info[tid0].bit & 1;
        for (int k = 1; k < it->second.size(); ++k)
          {
          if ((info[it->second[k]].bit & 1) == bit)
            ds.union_sets(tid0, it->second[k]);
          else
            {
            if (tid0_ == -1)
              tid0_ = int(it->second[k]);
            else
              ds.union_sets(tid0_, it->second[k]);
            }
          }
        }
      }

    // we re-organize the results of the union find as follows:
    std::vector<uint32_t> uq_ids(triangles.size(), (uint32_t)-1);
    std::vector< std::vector<uint32_t> > components;
    for (uint32_t i = 0; i < (uint32_t)triangles.size(); ++i)
      {
      auto ufid = ds.find(i);
      if (uq_ids[ufid] == (uint32_t)-1)
        { // unassigned
        uint32_t N = (uint32_t)components.size();
        components.push_back(std::vector<uint32_t>());

        uq_ids[ufid] = uq_ids[i] = N;
        components[N].push_back(i);
        }
      else { // assigned already
        uq_ids[i] = uq_ids[ufid]; // propagate assignment
        components[uq_ids[i]].push_back(i);
        }
      }
    return components;

    }

  bool tria_intersects_z_plane_through_point(std::vector<jtk::vec2<uint32_t>>& intersecting_edges, const jtk::vec3<double>& pt, const jtk::vec3<uint32_t>& tria, const std::vector<jtk::vec3<double>>& vert)
    {
    intersecting_edges.clear();
    for (uint32_t i = 0; i < 3; ++i)
      {
      jtk::vec2<uint32_t> e(tria[i], tria[(i + 1) % 3]);
      if (e[1] < e[0])
        std::swap(e[1], e[0]);
      auto V0 = vert[e[0]];
      auto V1 = vert[e[1]];
      if (V0[2] < pt[2])
        if (V1[2] > pt[2])
          intersecting_edges.push_back(e);
      if (V0[2] > pt[2])
        if (V1[2] < pt[2])
          intersecting_edges.push_back(e);
      if (V0[2] == pt[2] || V1[2] == pt[2])
        intersecting_edges.push_back(e);
      }
    return !intersecting_edges.empty();
    }


  bool _intersects(const jtk::vec3<double>& ray_origin, const jtk::vec3<double>& ray_direction, double ray_min, double ray_max, const jtk::vec3<double>& V0, const jtk::vec3<double>& V1, const jtk::vec3<double>& V2)
    { // Moeller-Trumbore algorithm
    auto edge1 = V1 - V0;
    auto edge2 = V2 - V0;
    auto pvec = cross(ray_direction, edge2);
    auto det = dot(edge1, pvec);
    if (det == 0.0)
      return false;
    double invDet = 1.0 / det;
    auto tvec = ray_origin - V0;
    std::array<double, 3> bary;
    bary[1] = dot(tvec, pvec) * invDet;
    if (bary[1] < 0 || bary[1] > 1)
      return false;
    auto qvec = cross(tvec, edge1);
    bary[2] = dot(ray_direction, qvec) * invDet;
    if (bary[2] < 0 || bary[1] + bary[2] > 1)
      return false;
    double t = dot(edge2, qvec) * invDet;
    //bary[0] = 1 - bary[1] - bary[2];
    if (t < ray_min)
      return false;
    if (t > ray_max)
      return false;
    return true;
    }

  static const double dirx[6] = { 1.0, 0.0, 0.0, -1.0, 0.0, 0.0 };
  static const double diry[6] = { 0.0, 1.0, 0.0, 0.0, -1.0, 0.0 };
  static const double dirz[6] = { 0.0, 0.0, 1.0, 0.0, 0.0, -1.0 };

  bool is_inside(const jtk::vec3<double>& pt, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<double>>& vertices, const std::vector<triangle_info>& info, const aabvh<uint32_t>& tria_aabvh, uint8_t mesh_id)
    {
    auto dir = normalize(jtk::vec3<double>(drand(0.5, 1.5), drand(0.5, 1.5), drand(0.5, 1.5)));
    int insides = 0;

    {

    int winding = 0;

    double ray_min = 0.0;
    double ray_max = std::numeric_limits<double>::infinity();

    tria_aabvh.for_each_through_ray(pt, dir, ray_min, ray_max, [&](uint32_t t)
      {
      if ((info[t].bit & 1) == mesh_id)
        return;
      auto tria = triangles[t];
      uint32_t a = tria[0];
      uint32_t b = tria[1];
      uint32_t c = tria[2];
      double flip = 1.0;
      // normalize vertex order (to prevent leaks)
      if (a > b)
        {
        std::swap(a, b);
        flip = -flip;
        }
      if (b > c)
        {
        std::swap(b, c);
        flip = -flip;
        }
      if (a > b)
        {
        std::swap(a, b);
        flip = -flip;
        }
      auto va = vertices[a];
      auto vb = vertices[b];
      auto vc = vertices[c];
      if (_intersects(pt, dir, ray_min, ray_max, va, vb, vc))
        {

        auto normal = flip * cross(vb - va, vc - va);
        if (dot(normal, dir) > 0.0)
          ++winding;
        else
          --winding;
        };
      });

    if (winding % 2 == 1)
      ++insides;

    }
    return insides > 0;
    }

  void _classify_triangles(std::vector<triangle_info>& info, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<double>>& vertices, const cork_options& options)
    {
    jtk::timer tim1;
    tim1.start();
    e_t_map edge_map = edges_to_triangles_map((uint32_t)vertices.size(), triangles);
    jtk::timer tim2;
    tim2.start();
    auto components = _compute_components(info, triangles, edge_map);
    if (options.p_str)
      {
      *options.p_str << "[LIBCORK] " << "Found " << (int)components.size() << " components\n";
      *options.p_str << "[LIBCORK] " << "Computing components took " << tim2.time_elapsed() << "s\n";
      }
    jtk::timer tim5;
    tim5.start();
    auto tria_aabvh = compute_triangle_aabvh(triangles, vertices);
    if (options.p_str)
      *options.p_str << "[LIBCORK] " << "Computing triangle aabvh took " << tim5.time_elapsed() << "s\n";
    std::vector<bool> visited(triangles.size(), false);
    _parallel_for((uint32_t)0, (uint32_t)components.size(), [&](uint32_t c)
      {
      const auto& comp = components[c];
      for (uint32_t id = 0; id < (uint32_t)comp.size(); ++id)
        {
        if (!visited[comp[id]])
          {
          uint8_t mesh_id = info[comp[id]].bit & 1;
          auto P = (vertices[triangles[comp[id]][0]] + vertices[triangles[comp[id]][1]] + vertices[triangles[comp[id]][2]]) / 3.0;
          bool inside = is_inside(P, triangles, vertices, info, tria_aabvh, mesh_id);
          info[comp[id]].bit |= inside ? 2 : 0;
          visited[comp[id]] = true;

          std::queue<uint32_t> qu;
          qu.push(comp[id]);
          while (!qu.empty())
            {
            auto t = qu.front();
            qu.pop();
            if (info[t].marked)
              continue;
            for (int i = 0; i < 3; ++i)
              {
              int j = (i + 1) % 3;
              jtk::vec2<uint32_t> e(triangles[t][i], triangles[t][j]);
              if (e[1] < e[0])
                std::swap(e[1], e[0]);
              int inside_sig = info[t].bit & 2;
              if (is_edge_with_marked_triangle(e, info, edge_map))
                continue;
              if (is_edge_on_intersection(e, info, edge_map))
                {
                inside_sig ^= 2;
                }

              for (auto t2 : edge_map[e[0]][e[1]])
                {
                if (visited[t2])
                  continue;
                if ((info[t2].bit & 1) != mesh_id)
                  {
                  continue;
                  }

                info[t2].bit |= inside_sig;

                visited[t2] = true;
                qu.push(t2);
                }
              }
            }
          }
        }

      }, options);
    if (options.p_str)
      *options.p_str << "[LIBCORK] " << "Classifying triangles took " << tim1.time_elapsed() << "s\n";
    }

  enum ETriangleOperation
    {
    DELETE_TRIANGLE,
    FLIP_TRIANGLE,
    KEEP_TRIANGLE
    };

  void _DeleteAndFlip(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<double>>& /*vertices*/, const std::vector<triangle_info>& info, std::function<ETriangleOperation(int bit)> classify)
    {
    std::vector<uint32_t> to_delete;
    std::vector<uint32_t> to_flip;
    for (uint32_t t = 0; t < (uint32_t)triangles.size(); ++t)
      {
      auto tria_op = classify(info[t].bit);
      switch (tria_op)
        {
        case DELETE_TRIANGLE: to_delete.push_back(t); break;
        case FLIP_TRIANGLE: to_flip.push_back(t); break;
        case KEEP_TRIANGLE: break;
        default: break;
        }
      }
    for (auto t : to_flip)
      {
      auto& tria = triangles[t];
      std::swap(tria[0], tria[2]);
      }
    jtk::delete_items(triangles, to_delete);
    }

  }

bool has_intersections(const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& vertices, const cork_options& options)
  {
  std::atomic<bool> degen{ false };
  std::atomic<bool> found_intersection{ false };
  auto verticesd = vertices_to_double(vertices);
  quantization q(get_magnitude(verticesd), 30);
  verticesd = get_scaled_vertices(verticesd, q);
  auto edge_to_tria_map = edges_to_triangles_map((uint32_t)vertices.size(), triangles);
  auto edge_aabvh = compute_edge_aabvh(edge_to_tria_map, verticesd);

  aabvh_edge_triangle_parallel([&](const jtk::vec2<uint32_t>& e, const jtk::vec3<uint32_t>& t, uint32_t)->bool
    {
    bool deg;
    if (intersects(deg, e, t, verticesd, q))
      {
      found_intersection = true;
      return false; // break from aabvh_edge_triangle
      }
    if (deg)
      {
      degen = true;
      return false;
      }
    return true; // continue with aabvh_edge_triangle
    }, edge_aabvh, triangles, verticesd, options);
  if (degen)
    {
    if (options.p_str)
      *options.p_str << "[LIBCORK] " << "Degeneracies were detected.\n";
    }
  return found_intersection;
  }

bool is_closed(const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& vertices)
  {
  std::vector<std::map<uint32_t, uint32_t > > edge_count(vertices.size());
  for (uint32_t t = 0; t < (uint32_t)triangles.size(); ++t)
    {
    for (uint32_t i = 0; i < 3; ++i)
      {
      uint32_t j = (i + 1) % 3;
      uint32_t v0 = triangles[t][i];
      uint32_t v1 = triangles[t][j];
      auto find_edge = edge_count[v0].find(v1);
      if (find_edge == edge_count[v0].end())
        edge_count[v0][v1] = 1;
      else
        ++edge_count[v0][v1];
      find_edge = edge_count[v1].find(v0);
      if (find_edge == edge_count[v1].end())
        edge_count[v1][v0] = 0;
      }
    }
  for (uint32_t i = 0; i < (uint32_t)edge_count.size(); ++i)
    {
    for (auto it = edge_count[i].begin(); it != edge_count[i].end(); ++it)
      {
      if (it->second != edge_count[it->first][i])
        return false;
      }
    }
  return true;
  }

bool is_solid(const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& vertices, const cork_options& options)
  {
  if (!is_closed(triangles, vertices))
    return false;
  if (has_intersections(triangles, vertices, options))
    return false;
  return true;
  }

void compute_union(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles_left, const std::vector<jtk::vec3<float>>& vertices_left, const float* cs_left, const std::vector<jtk::vec3<uint32_t>>& triangles_right, const std::vector<jtk::vec3<float>>& vertices_right, const float* cs_right, const cork_options& options)
  {
  std::vector<triangle_info> info;
  std::vector<jtk::vec3<double>> verticesd;
  _resolve_intersections(verticesd, triangles, info, triangles_left, vertices_left, cs_left, triangles_right, vertices_right, cs_right, options);
  _classify_triangles(info, triangles, verticesd, options);
  _DeleteAndFlip(triangles, verticesd, info, [](int bit) -> ETriangleOperation
    {
    if ((bit & 2) == 2)
      return DELETE_TRIANGLE;
    else
      return KEEP_TRIANGLE;
    });
  vertices = vertices_to_float(verticesd);
  }

void compute_difference(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles_left, const std::vector<jtk::vec3<float>>& vertices_left, const float* cs_left, const std::vector<jtk::vec3<uint32_t>>& triangles_right, const std::vector<jtk::vec3<float>>& vertices_right, const float* cs_right, const cork_options& options)
  {
  std::vector<triangle_info> info;
  std::vector<jtk::vec3<double>> verticesd;
  _resolve_intersections(verticesd, triangles, info, triangles_left, vertices_left, cs_left, triangles_right, vertices_right, cs_right, options);
  _classify_triangles(info, triangles, verticesd, options);
  _DeleteAndFlip(triangles, verticesd, info, [](int bit) -> ETriangleOperation
    {
    if (bit == 2 || bit == 1)
      return DELETE_TRIANGLE;
    else if (bit == 3)
      return FLIP_TRIANGLE;
    else
      return KEEP_TRIANGLE;
    });
  vertices = vertices_to_float(verticesd);
  }


void compute_intersection(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles_left, const std::vector<jtk::vec3<float>>& vertices_left, const float* cs_left, const std::vector<jtk::vec3<uint32_t>>& triangles_right, const std::vector<jtk::vec3<float>>& vertices_right, const float* cs_right, const cork_options& options)
  {
  std::vector<triangle_info> info;
  std::vector<jtk::vec3<double>> verticesd;
  _resolve_intersections(verticesd, triangles, info, triangles_left, vertices_left, cs_left, triangles_right, vertices_right, cs_right, options);
  _classify_triangles(info, triangles, verticesd, options);
  _DeleteAndFlip(triangles, verticesd, info, [](int bit) -> ETriangleOperation
    {
    if ((bit & 2) == 0)
      return DELETE_TRIANGLE;
    else
      return KEEP_TRIANGLE;
    });
  vertices = vertices_to_float(verticesd);
  }


void compute_symmetric_difference(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles_left, const std::vector<jtk::vec3<float>>& vertices_left, const float* cs_left, const std::vector<jtk::vec3<uint32_t>>& triangles_right, const std::vector<jtk::vec3<float>>& vertices_right, const float* cs_right, const cork_options& options)
  {
  std::vector<triangle_info> info;
  std::vector<jtk::vec3<double>> verticesd;
  _resolve_intersections(verticesd, triangles, info, triangles_left, vertices_left, cs_left, triangles_right, vertices_right, cs_right, options);
  _classify_triangles(info, triangles, verticesd, options);
  _DeleteAndFlip(triangles, verticesd, info, [](int bit) -> ETriangleOperation
    {
    if ((bit & 2) == 0)
      return KEEP_TRIANGLE;
    else
      return FLIP_TRIANGLE;
    });
  vertices = vertices_to_float(verticesd);
  }

void resolve_intersections(std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles_in, const std::vector<jtk::vec3<float>>& vertices_in, const cork_options& options)
  {
  std::vector<triangle_info> info;
  std::vector<jtk::vec3<double>> verticesd = vertices_to_double(vertices_in);
  triangles = triangles_in;
  _resolve_intersections(info, triangles, verticesd, true, options);
  vertices = vertices_to_float(verticesd);
  }


diagnostics diagnose(const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<float>>& vertices, const cork_options& options)
  {
  diagnostics d;
  std::vector<jtk::vec3<double>> verticesd = vertices_to_double(vertices);

  jtk::vec3<double> center(0, 0, 0);
  if (options.center)
    {
    for (const auto& pt : verticesd)
      center = center + pt;
    center = center / (double)verticesd.size();
    for (auto& pt : verticesd)
      {
      pt = pt - center;
      }
    }

  quantization q(get_magnitude(verticesd), 30);
  verticesd = get_scaled_vertices(verticesd, q);

  d.number_of_vertices = (uint32_t)vertices.size();
  d.number_of_triangles = (uint32_t)triangles.size();

  std::atomic<int> degen = 0;
  std::atomic<int> isct = 0;
  auto edge_to_tria_map = edges_to_triangles_map((uint32_t)vertices.size(), triangles);
  auto edge_aabvh = compute_edge_aabvh(edge_to_tria_map, verticesd);
  aabvh_edge_triangle_parallel_allow_degeneracies([&](const jtk::vec2<uint32_t>& e, const jtk::vec3<uint32_t>& t, uint32_t /*tria_index*/)->bool
    {
    bool deg = false;
    if (intersects(deg, e, t, verticesd, q))
      {
      ++isct;
      }
    if (deg)
      {
      ++degen;
      }
    return true;
    }, edge_aabvh, triangles, verticesd, options);

  d.intersections = (uint32_t)isct;
  d.degenerate = (uint32_t)degen;
  return d;
  }

