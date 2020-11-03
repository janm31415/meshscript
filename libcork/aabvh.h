#pragma once

#include <vector>
#include <functional>
#include <algorithm>

#include <jtk/vec.h>

#include <cassert>

/*
Axis Aligned Bounding Volume Hierarchy
*/

#undef min
#undef max

inline jtk::boundingbox3d<double> make_boundingbox(const jtk::vec3<double>& _min, const jtk::vec3<double>& _max)
  {
  jtk::boundingbox3d<double> bb;
  bb.min = jtk::min(_min, _max);
  bb.max = jtk::max(_min, _max);
  return bb;
  }

inline jtk::boundingbox3d<double> make_boundingbox(double x, double y, double z, double X, double Y, double Z)
  {
  jtk::boundingbox3d<double> bb;
  bb.min = jtk::min(jtk::vec3<double>(x, y, z), jtk::vec3<double>(X, Y, Z));
  bb.max = jtk::max(jtk::vec3<double>(x, y, z), jtk::vec3<double>(X, Y, Z));
  return bb;
  }

template <class GeomIdx>
struct aabvh_entity
  {
  jtk::boundingbox3d<double> bbox;
  jtk::vec3<double> point;
  GeomIdx idx;
  };

template <class GeomIdx>
struct aabb_node
  {
  jtk::boundingbox3d<double> bbox;
  int left;
  int right;
  std::vector<GeomIdx> aabvh_idxs;
  std::vector<jtk::boundingbox3d<double>> aabvh_bboxs;
  inline bool is_leaf() const
    {
    return left < 0;
    }
  };

inline bool intersects(const jtk::boundingbox3d<double>& left, const jtk::boundingbox3d<double>& right)
  {
  return left.min[0] <= right.max[0] && left.max[0] >= right.min[0] &&
    left.min[1] <= right.max[1] && left.max[1] >= right.min[1] &&
    left.min[2] <= right.max[2] && left.max[2] >= right.min[2];
  }

inline bool intersects(const jtk::boundingbox3d<double>& bb, const jtk::vec3<double>& origin, const jtk::vec3<double>& direction, double tmin, double tmax)
  {
  double txmin = (bb.min[0] - origin[0]) / direction[0]; // Division by zero is allowed
  double txmax = (bb.max[0] - origin[0]) / direction[0]; // Division by zero is allowed
  if (txmin > txmax)
    std::swap(txmin, txmax);
  double tymin = (bb.min[1] - origin[1]) / direction[1]; // Division by zero is allowed
  double tymax = (bb.max[1] - origin[1]) / direction[1]; // Division by zero is allowed
  if (tymin > tymax)
    std::swap(tymin, tymax);
  if ((txmin > tymax) || (tymin > txmax))
    return false;
  if (tymin > txmin)
    txmin = tymin;
  if (tymax < txmax)
    txmax = tymax;
  double tzmin = (bb.min[2] - origin[2]) / direction[2]; // Division by zero is allowed
  double tzmax = (bb.max[2] - origin[2]) / direction[2]; // Division by zero is allowed
  if (tzmin > tzmax)
    std::swap(tzmin, tzmax);
  if ((txmin > tzmax) || (tzmin > txmax))
    return false;
  if (tzmin > txmin)
    txmin = tzmin;
  if (tzmax < txmax)
    txmax = tzmax;
  if (txmin > tmax)
    return false;
  if (txmax < tmin)
    return false;
  return true;
  }

template <class GeomIdx, int leaf_size = 8>
class aabvh
  {
  typedef typename std::vector<aabvh_entity<GeomIdx>>::iterator iterator;
  public:
    aabvh(std::vector<aabvh_entity<GeomIdx>>& geometries)
      {
      assert(!geometries.empty());
      root_id = construct_tree(geometries.begin(), geometries.end(), 0);
      }

    ~aabvh() {}

    void for_each_in_box(const jtk::boundingbox3d<double>& bbox,
      std::function< void(const GeomIdx& idx) > action) const
      {
      std::vector<size_t> stack;
      stack.push_back(root_id);
      while (!stack.empty())
        {
        const auto& node = nodes[stack.back()];
        stack.pop_back();
        if (!intersects(node.bbox, bbox))
          continue;
        if (node.is_leaf())
          {
          for (int i = 0; i < node.aabvh_idxs.size(); ++i)
            {
            if (intersects(node.aabvh_bboxs[i], bbox))
              action(node.aabvh_idxs[i]);
            }
          }
        else
          {
          stack.push_back(node.left);
          stack.push_back(node.right);
          }
        }
      }

    void for_each_through_ray(const jtk::vec3<double>& ray_origin, const jtk::vec3<double>& ray_direction,
      double ray_min, double ray_max, std::function< void(const GeomIdx& idx) > action) const
      {
      std::vector<size_t> stack;
      stack.push_back(root_id);
      while (!stack.empty())
        {
        const auto& node = nodes[stack.back()];
        stack.pop_back();
        if (!intersects(node.bbox, ray_origin, ray_direction, ray_min, ray_max))
          continue;
        if (node.is_leaf())
          {
          for (int i = 0; i < node.aabvh_idxs.size(); ++i)
            {
            if (intersects(node.aabvh_bboxs[i], ray_origin, ray_direction, ray_min, ray_max))
              action(node.aabvh_idxs[i]);
            }
          }
        else
          {
          stack.push_back(node.left);
          stack.push_back(node.right);
          }
        }
      }

    void for_each_through_ray(const jtk::vec3<float>& ray_origin, const jtk::vec3<float>& ray_direction,
      float ray_min, float ray_max, std::function< void(const GeomIdx& idx) > action) const
      {
      std::vector<size_t> stack;
      stack.push_back(root_id);
      while (!stack.empty())
        {
        const auto& node = nodes[stack.back()];
        stack.pop_back();
        if (!intersects(node.bbox, ray_origin, ray_direction, ray_min, ray_max))
          continue;
        if (node.is_leaf())
          {
          for (int i = 0; i < node.aabvh_idxs.size(); ++i)
            {
            if (intersects(node.aabvh_bboxs[i], ray_origin, ray_direction, ray_min, ray_max))
              action(node.aabvh_idxs[i]);
            }
          }
        else
          {
          stack.push_back(node.left);
          stack.push_back(node.right);
          }
        }
      }

    size_t number_of_leafs() const
      {
      return leaf_positions.size();
      }

    const std::vector<GeomIdx>& get_leaf_entities(size_t leaf_idx) const
      {
      return nodes[leaf_positions[leaf_idx]].aabvh_idxs;
      }

    const jtk::boundingbox3d<double>& get_leaf_bbox(size_t leaf_idx) const
      {
      return nodes[leaf_positions[leaf_idx]].bbox;
      }

  private:
    int construct_tree(iterator first, iterator last, int dim)
      {
      assert(first != last);
      auto sz = std::distance(first, last);
      if (sz <= leaf_size)
        {
        aabb_node<GeomIdx> n;
        n.left = -1;
        n.aabvh_idxs.resize(sz);
        n.aabvh_bboxs.resize(sz);
        n.bbox = first->bbox;
        auto it = first;
        for (int i = 0; i < int(sz); ++i, ++it)
          {
          n.aabvh_idxs[i] = it->idx;
          n.aabvh_bboxs[i] = it->bbox;
          n.bbox = unite(n.bbox, it->bbox);
          }
        int node_id = int(nodes.size());
        nodes.push_back(n);
        leaf_positions.push_back(node_id);
        return node_id;
        }
      iterator mid = first + (last - first) / 2;
      std::nth_element(first, mid, last, [&](const aabvh_entity<GeomIdx>& left, const aabvh_entity<GeomIdx>& right)
        {
        return left.point[dim] < right.point[dim];
        });
      aabb_node<GeomIdx> n;
      assert(first != mid);
      n.left = construct_tree(first, mid, (dim + 1) % 3);
      n.right = construct_tree(mid, last, (dim + 1) % 3);
      n.bbox = unite(nodes[n.left].bbox, nodes[n.right].bbox);
      int node_id = int(nodes.size());
      nodes.push_back(n);
      return node_id;
      }

  private:
    int root_id;
    std::vector<aabb_node<GeomIdx>> nodes;
    std::vector<size_t> leaf_positions;
  };
