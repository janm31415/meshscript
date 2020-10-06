#pragma once

#include <jtk/vec.h>
#include <vector>
#include <limits>
#include <set>
#include <algorithm>

namespace jtk
  {

  namespace point_tree_details
    {

    template <class _Ty>
    inline _Ty _Sqr(_Ty _Value)
      {
      return _Value * _Value;
      }

    template <class _Ty, class _Point_type, int _Dim>
    struct _Dist_sqr_impl
      {
      static _Ty fun(const _Point_type& _Pt_1, const _Point_type& _Pt_2)
        {
        _Ty dist = _Sqr(_Pt_1[0] - _Pt_2[0]);
        for (size_t i = 1; i < _Dim; ++i)
          dist += _Sqr(_Pt_1[i] - _Pt_2[i]);
        return dist;
        }
      };

    template <class _Ty, class _Point_type>
    struct _Dist_sqr_impl<_Ty, _Point_type, 1>
      {
      static _Ty fun(const _Point_type& _Pt_1, const _Point_type& _Pt_2)
        {
        return _Sqr(_Pt_1[0] - _Pt_2[0]);
        }
      };

    template <class _Ty, class _Point_type>
    struct _Dist_sqr_impl<_Ty, _Point_type, 2>
      {
      static _Ty fun(const _Point_type& _Pt_1, const _Point_type& _Pt_2)
        {
        return _Sqr(_Pt_1[0] - _Pt_2[0]) + _Sqr(_Pt_1[1] - _Pt_2[1]);
        }
      };

    template <class _Ty, class _Point_type>
    struct _Dist_sqr_impl<_Ty, _Point_type, 3>
      {
      static _Ty fun(const _Point_type& _Pt_1, const _Point_type& _Pt_2)
        {
        return _Sqr(_Pt_1[0] - _Pt_2[0]) + _Sqr(_Pt_1[1] - _Pt_2[1]) + _Sqr(_Pt_1[2] - _Pt_2[2]);
        }
      };

    template <class _Ty, class _Point_type>
    struct _Dist_sqr_impl<_Ty, _Point_type, 4>
      {
      static _Ty fun(const _Point_type& _Pt_1, const _Point_type& _Pt_2)
        {
        return _Sqr(_Pt_1[0] - _Pt_2[0]) + _Sqr(_Pt_1[1] - _Pt_2[1]) + _Sqr(_Pt_1[2] - _Pt_2[2]) + _Sqr(_Pt_1[3] - _Pt_2[3]);
        }
      };

    template <class _Ty, class _Point_type, int _Dim>
    inline _Ty _Dist_sqr(const _Point_type& _Pt_1, const _Point_type& _Pt_2)
      {
      return _Dist_sqr_impl<_Ty, _Point_type, _Dim>::fun(_Pt_1, _Pt_2);
      }

    template <class _Ty, class _Point_type, int _Dim>
    struct aabb_node
      {
      _Point_type _Min, _Max;

      aabb_node() {}

      aabb_node(const _Point_type& _P0) : _Min(_P0), _Max(_P0)
        {
        }

      aabb_node(const _Point_type& _V0, const _Point_type& _V1, const _Point_type& _V2) : _Min(_V0), _Max(_V1)
        {
        for (size_t i = 0; i < _Dim; ++i)
          {
          if (_Min[i] > _Max[i])
            std::swap(_Min[i], _Max[i]);
          }
        add_point(_V2);
        }

      aabb_node(const aabb_node& _Left, const aabb_node& _Right) : _Min(_Left._Min), _Max(_Left._Max)
        {
        add_point(_Right._Min);
        add_point(_Right._Max);
        }

      void add_point(const _Point_type& _Pt)
        {
        for (size_t i = 0; i < _Dim; ++i)
          {
          _Min[i] = std::min<_Ty>(_Min[i], _Pt[i]);
          _Max[i] = std::max<_Ty>(_Max[i], _Pt[i]);
          }
        }

      _Point_type& Border(size_t i)
        {
        return (&_Min)[i];
        }

      const _Point_type& Border(size_t i) const
        {
        return (&_Min)[i];
        }

      template <int _D>
      inline _Ty distance_sqr(const _Point_type& _Pt) const
        {
        _Ty dist = _Sqr(_Pt[0] - std::min<_Ty>(_Max[0], std::max<_Ty>(_Min[0], _Pt[0])));
        for (size_t i = 1; i < _D; ++i)
          dist += _Sqr(_Pt[i] - std::min<_Ty>(_Max[i], std::max<_Ty>(_Min[i], _Pt[i])));
        return dist;
        }

      template <>
      inline _Ty distance_sqr<4>(const _Point_type& _Pt) const
        {
        return _Sqr(_Pt[0] - std::min<_Ty>(_Max[0], std::max<_Ty>(_Min[0], _Pt[0])))
          + _Sqr(_Pt[1] - std::min<_Ty>(_Max[1], std::max<_Ty>(_Min[1], _Pt[1])))
          + _Sqr(_Pt[2] - std::min<_Ty>(_Max[2], std::max<_Ty>(_Min[2], _Pt[2])))
          + _Sqr(_Pt[3] - std::min<_Ty>(_Max[3], std::max<_Ty>(_Min[3], _Pt[3])));
        }

      template <>
      inline _Ty distance_sqr<3>(const _Point_type& _Pt) const
        {
        return _Sqr(_Pt[0] - std::min<_Ty>(_Max[0], std::max<_Ty>(_Min[0], _Pt[0])))
          + _Sqr(_Pt[1] - std::min<_Ty>(_Max[1], std::max<_Ty>(_Min[1], _Pt[1])))
          + _Sqr(_Pt[2] - std::min<_Ty>(_Max[2], std::max<_Ty>(_Min[2], _Pt[2])));
        }

      template <>
      inline _Ty distance_sqr<2>(const _Point_type& _Pt) const
        {
        return _Sqr(_Pt[0] - std::min<_Ty>(_Max[0], std::max<_Ty>(_Min[0], _Pt[0])))
          + _Sqr(_Pt[1] - std::min<_Ty>(_Max[1], std::max<_Ty>(_Min[1], _Pt[1])));
        }

      template <>
      inline _Ty distance_sqr<1>(const _Point_type& _Pt) const
        {
        return _Sqr(_Pt[0] - std::min<_Ty>(_Max[0], std::max<_Ty>(_Min[0], _Pt[0])));
        }

      inline _Ty center(size_t dim) const
        {
        return (_Min[dim] + _Max[dim]) / 2.0;
        }

      };

    template <class _Ty, class _Point_type, int _Dim>
    struct kd_point_node
      {
      aabb_node<_Ty, _Point_type, _Dim> _Node;
      size_t _Child_1, _Child_2;
      };

    template <class _Point_type>
    struct point_compare
      {
      public:
        point_compare(size_t dim)
          : _Dim(dim) {}

        inline bool operator()(const _Point_type& _Left, const _Point_type& _Right) const
          {
          return _Left[_Dim] < _Right[_Dim];
          }

      private:
        size_t _Dim;
      };

    template<class _Ty, class _Result_type>
    struct k_set_compare
      {
      bool operator () (const std::pair<_Ty, _Result_type>& _Left, const std::pair<_Ty, _Result_type>& _Right)
        {
        return _Left.first < _Right.first;
        }
      };

    }

  struct default_point_tree_traits
    {
    typedef double value_type;
    enum { dimension = 3 };
    typedef vec3<double> point_type;
    };

  template <class _Traits = default_point_tree_traits>
  class point_tree
    {
    public:
      typedef typename _Traits::point_type point_type;
      typedef typename _Traits::value_type value_type;

      point_tree();

      typename _Traits::point_type find_nearest(value_type& distance, const typename _Traits::point_type& point) const;
      template <class Predicate>
      typename _Traits::point_type find_nearest_if(value_type& distance, const typename _Traits::point_type& point, Predicate pred) const;

      std::vector<typename _Traits::point_type> find_nearest_within_radius(const value_type& radius, const typename _Traits::point_type& point) const;

      std::vector<typename _Traits::point_type> find_k_nearest(int k, const typename _Traits::point_type& point) const;

      template <class Predicate>
      std::vector<typename _Traits::point_type> find_k_nearest_if(int k, const typename _Traits::point_type& point, Predicate pred) const;

      template <class Iter>
      void build_tree(Iter _First, Iter _Last);

      // efficient_build_tree will rearrange the order of points that are given as input
      template <class Iter>
      void efficient_build_tree(Iter _First, Iter _Last);

      bool empty() const;

    private:

      template <class Iter>
      void _Init(Iter _First, Iter _Last);

      template <class Iter>
      size_t _Optimise(Iter _First, Iter _Last, size_t _Level = 0);

    private:
      size_t _Size;
      std::vector<typename point_tree_details::kd_point_node<value_type, point_type, _Traits::dimension>> _Kd;
    };

  template <class _Traits>
  point_tree<_Traits>::point_tree() :
    _Size(0)
    {
    }

  template <class _Traits>
  template <class Iter>
  size_t point_tree<_Traits>::_Optimise(Iter _First, Iter _Last, size_t _Level)
    {
    using namespace point_tree_details;
    kd_point_node<value_type, point_type, _Traits::dimension> _N;
    point_compare<point_type> _Compare(_Level % _Traits::dimension);
    Iter _Mid = _First + (_Last - _First) / 2;
    std::nth_element(_First, _Mid, _Last, _Compare);
    if (_Mid != _First)
      {
      _N._Child_1 = _Optimise(_First, _Mid, _Level + 1);
      _N._Child_2 = _Optimise(_Mid, _Last, _Level + 1);
      _N._Node = aabb_node<value_type, point_type, _Traits::dimension>(_Kd[_N._Child_1]._Node, _Kd[_N._Child_2]._Node);
      }
    else
      {
      _N._Child_1 = _N._Child_2 = size_t(-1);
      _N._Node = aabb_node<value_type, point_type, _Traits::dimension>(*_First);
      }
    _Kd.push_back(_N);
    return _Kd.size() - 1;
    }

  template <class _Traits>
  template <class Iter>
  void point_tree<_Traits>::_Init(Iter _First, Iter _Last)
    {
    using namespace point_tree_details;
    _Kd.clear();
    _Kd.reserve(_Size * 2);
    _Optimise(_First, _Last, 0);
    }

  template <class _Traits>
  bool point_tree<_Traits>::empty() const
    {
    return _Size == 0;
    }

  template <class _Traits>
  template <class Iter>
  void point_tree<_Traits>::build_tree(Iter _First, Iter _Last)
    {
    if (_First == _Last)
      return;
    std::vector<typename _Traits::point_type> _Leafs;
    _Size = std::distance(_First, _Last);
    _Leafs.resize(_Size);
    std::copy(_First, _Last, _Leafs.begin());
    _Init(_Leafs.begin(), _Leafs.end());
    }

  template <class _Traits>
  template <class Iter>
  void point_tree<_Traits>::efficient_build_tree(Iter _First, Iter _Last)
    {
    _Size = std::distance(_First, _Last);
    _Init(_First, _Last);
    }

  template <class _Traits>
  typename _Traits::point_type point_tree<_Traits>::find_nearest(value_type& distance, const typename _Traits::point_type& point) const
    {
    std::vector<size_t> _St;
    _St.reserve(size_t(std::log(_Size)*2.0));
    _St.push_back(_Kd.size() - 1);
    typename _Traits::point_type _Pt = _Kd[_St.back()]._Node._Min;
    value_type _Best_d = std::numeric_limits<value_type>::infinity();
    while (!_St.empty())
      {
      const auto& _Kd_node = _Kd[_St.back()];
      _St.pop_back();
      if (_Kd_node._Child_1 == size_t(-1)) // a leaf
        {
        value_type _Dist = point_tree_details::_Dist_sqr<value_type, point_type, _Traits::dimension>(point, _Kd_node._Node._Min);
        if (_Dist < _Best_d)
          {
          _Best_d = _Dist;
          _Pt = _Kd_node._Node._Min;
          }
        }
      else
        {
        value_type _Min_d1 = _Kd[_Kd_node._Child_1]._Node.distance_sqr<_Traits::dimension>(point);
        value_type _Min_d2 = _Kd[_Kd_node._Child_2]._Node.distance_sqr<_Traits::dimension>(point);

        if (_Min_d1 < _Min_d2)
          {
          if (_Best_d > _Min_d2)
            {
            _St.push_back(_Kd_node._Child_2);
            _St.push_back(_Kd_node._Child_1);
            }
          else if (_Best_d > _Min_d1)
            {
            _St.push_back(_Kd_node._Child_1);
            }
          }
        else
          {
          if (_Best_d > _Min_d1)
            {
            _St.push_back(_Kd_node._Child_1);
            _St.push_back(_Kd_node._Child_2);
            }
          else if (_Best_d > _Min_d2)
            {
            _St.push_back(_Kd_node._Child_2);
            }
          }
        }
      }
    distance = std::sqrt(_Best_d);
    return _Pt;
    }

  template <class _Traits>
  template <class Predicate>
  typename _Traits::point_type point_tree<_Traits>::find_nearest_if(value_type& distance, const typename _Traits::point_type& point, Predicate pred) const
    {
    typename _Traits::point_type _Pt;
    std::vector<size_t> _St;
    _St.reserve(size_t(std::log(_Size)*2.0));
    _St.push_back(_Kd.size() - 1);
    value_type _Best_d = std::numeric_limits<value_type>::infinity();
    while (!_St.empty())
      {
      const auto& _Kd_node = _Kd[_St.back()];
      _St.pop_back();
      if (_Kd_node._Child_1 == size_t(-1)) // a leaf
        {
        if (pred(_Kd_node._Node._Min))
          {
          value_type _Dist = point_tree_details::_Dist_sqr<value_type, point_type, _Traits::dimension>(point, _Kd_node._Node._Min);
          if (_Dist < _Best_d)
            {
            _Best_d = _Dist;
            _Pt = _Kd_node._Node._Min;
            }
          }
        }
      else
        {
        value_type _Min_d1 = _Kd[_Kd_node._Child_1]._Node.distance_sqr<_Traits::dimension>(point);
        value_type _Min_d2 = _Kd[_Kd_node._Child_2]._Node.distance_sqr<_Traits::dimension>(point);

        if (_Min_d1 < _Min_d2)
          {
          if (_Best_d > _Min_d2)
            {
            _St.push_back(_Kd_node._Child_2);
            _St.push_back(_Kd_node._Child_1);
            }
          else if (_Best_d > _Min_d1)
            {
            _St.push_back(_Kd_node._Child_1);
            }
          }
        else
          {
          if (_Best_d > _Min_d1)
            {
            _St.push_back(_Kd_node._Child_1);
            _St.push_back(_Kd_node._Child_2);
            }
          else if (_Best_d > _Min_d2)
            {
            _St.push_back(_Kd_node._Child_2);
            }
          }
        }
      }
    distance = std::sqrt(_Best_d);
    return _Pt;
    }

  template <class _Traits>
  std::vector<typename _Traits::point_type> point_tree<_Traits>::find_nearest_within_radius(const value_type& radius, const typename _Traits::point_type& point) const
    {
    std::vector<typename _Traits::point_type> _Pts;
    std::vector<size_t> _St;
    _St.reserve(size_t(std::log(_Size)*2.0));
    _St.push_back(_Kd.size() - 1);
    value_type _Rad_2 = radius * radius;
    while (!_St.empty())
      {
      const auto& _Kd_node = _Kd[_St.back()];
      _St.pop_back();
      if (_Kd_node._Child_1 == size_t(-1)) // a leaf
        {
        value_type _Dist = point_tree_details::_Dist_sqr<value_type, point_type, _Traits::dimension>(point, _Kd_node._Node._Min);
        if (_Dist <= _Rad_2)
          {
          _Pts.push_back(_Kd_node._Node._Min);
          }
        }
      else
        {
        value_type _Min_d1 = _Kd[_Kd_node._Child_1]._Node.distance_sqr<_Traits::dimension>(point);
        value_type _Min_d2 = _Kd[_Kd_node._Child_2]._Node.distance_sqr<_Traits::dimension>(point);
        if (_Rad_2 > _Min_d2)
          _St.push_back(_Kd_node._Child_2);
        if (_Rad_2 > _Min_d1)
          _St.push_back(_Kd_node._Child_1);
        }
      }
    return _Pts;
    }

  template <class _Traits>
  std::vector<typename _Traits::point_type> point_tree<_Traits>::find_k_nearest(int k, const typename _Traits::point_type& point) const
    {
    std::vector<typename _Traits::point_type> _Pts;
    std::multiset<std::pair<typename _Traits::value_type, typename _Traits::point_type>, point_tree_details::k_set_compare<typename _Traits::value_type, typename _Traits::point_type>> _Candidates;
    std::vector<size_t> _St;
    _St.reserve(size_t(std::log(_Size)*2.0));
    _St.push_back(_Kd.size() - 1);
    value_type _K_best_d = std::numeric_limits<value_type>::infinity();
    while (!_St.empty())
      {
      const auto& _Kd_node = _Kd[_St.back()];
      _St.pop_back();
      if (_Kd_node._Child_1 == size_t(-1)) // a leaf
        {
        value_type _Dist = point_tree_details::_Dist_sqr<value_type, point_type, _Traits::dimension>(point, _Kd_node._Node._Min);
        if (_Dist <= _K_best_d)
          {
          if (_Candidates.size() >= k)
            _Candidates.erase(--_Candidates.end());
          _Candidates.insert(std::pair<typename _Traits::value_type, typename _Traits::point_type>(_Dist, _Kd_node._Node._Min));
          _K_best_d = _Candidates.rbegin()->first;
          }
        }
      else
        {
        value_type _Min_d1 = _Kd[_Kd_node._Child_1]._Node.distance_sqr<_Traits::dimension>(point);
        value_type _Min_d2 = _Kd[_Kd_node._Child_2]._Node.distance_sqr<_Traits::dimension>(point);
        if (_K_best_d > _Min_d2 || _Candidates.size() < k)
          _St.push_back(_Kd_node._Child_2);
        if (_K_best_d > _Min_d1 || _Candidates.size() < k)
          _St.push_back(_Kd_node._Child_1);
        }
      }
    for (const auto& _Cand : _Candidates)
      _Pts.push_back(_Cand.second);
    return _Pts;
    }

  template <class _Traits>
  template <class Predicate>
  std::vector<typename _Traits::point_type> point_tree<_Traits>::find_k_nearest_if(int k, const typename _Traits::point_type& point, Predicate pred) const
    {
    std::vector<typename _Traits::point_type> _Pts;
    std::multiset<std::pair<typename _Traits::value_type, typename _Traits::point_type>, details::k_set_compare<typename _Traits::value_type, typename _Traits::point_type>> _Candidates;
    std::vector<size_t> _St;
    _St.reserve(size_t(std::log(_Size)*2.0));
    _St.push_back(_Kd.size() - 1);
    value_type _K_best_d = std::numeric_limits<value_type>::infinity();
    while (!_St.empty())
      {
      const auto& _Kd_node = _Kd[_St.back()];
      _St.pop_back();
      if (_Kd_node._Child_1 == size_t(-1)) // a leaf
        {
        if (pred(_Kd_node._Node._Min))
          {
          value_type _Dist = point_tree_details::_Dist_sqr<value_type, point_type, _Traits::dimension>(point, _Kd_node._Node._Min);
          if (_Dist <= _K_best_d)
            {
            if (_Candidates.size() >= k)
              _Candidates.erase(--_Candidates.end());
            _Candidates.insert(std::pair<typename _Traits::value_type, typename _Traits::point_type>(_Dist, _Kd_node._Node._Min));
            _K_best_d = _Candidates.rbegin()->first;
            }
          }
        }
      else
        {
        value_type _Min_d1 = _Kd[_Kd_node._Child_1]._Node.distance_sqr<_Traits::dimension>(point);
        value_type _Min_d2 = _Kd[_Kd_node._Child_2]._Node.distance_sqr<_Traits::dimension>(point);
        if (_K_best_d > _Min_d2 || _Candidates.size() < k)
          _St.push_back(_Kd_node._Child_2);
        if (_K_best_d > _Min_d1 || _Candidates.size() < k)
          _St.push_back(_Kd_node._Child_1);
        }
      }
    for (const auto& _Cand : _Candidates)
      _Pts.push_back(_Cand.second);
    return _Pts;
    }

  }