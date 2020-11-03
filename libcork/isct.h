#pragma once

#include "ga4.h"

#include <cmath>
#include <limits>

namespace ga
  {
  /*
  The join of blades A and B is the span of their subspaces. if A /cap B is empty, then their join
  equals the outer product A ^ B.
  If A /cap B is not empty, then the following join functions will return 0, which is not necessarly the correct join, but for
  the purpose of the intersection algorithms this is sufficient.
  */

  template <class TType>
  bivector<TType> join(const vector<TType>& v1, const vector<TType>& v2)
    {
    return outer(v1, v2);
    }

  template <class TType>
  trivector<TType> join(const bivector<TType>& v1, const vector<TType>& v2)
    {
    return outer(v1, v2);
    }

  template <class TType>
  trivector<TType> join(const vector<TType>& v1, const bivector<TType>& v2)
    {
    return outer(v1, v2);
    }

  template <class TType>
  bivector<TType> join_abs(const vector<TType>& v1, const vector<TType>& v2)
    {
    bivector<TType> res;
    res.e0e1 = fabs(v1.e0 * v2.e1) + fabs(v2.e0 * v1.e1);
    res.e0e2 = fabs(v1.e0 * v2.e2) + fabs(v2.e0 * v1.e2);
    res.e0e3 = fabs(v1.e0 * v2.e3) + fabs(v2.e0 * v1.e3);
    res.e1e2 = fabs(v1.e1 * v2.e2) + fabs(v2.e1 * v1.e2);
    res.e1e3 = fabs(v1.e1 * v2.e3) + fabs(v2.e1 * v1.e3);
    res.e2e3 = fabs(v1.e2 * v2.e3) + fabs(v2.e2 * v1.e3);
    return res;
    }

  template <class TType>
  trivector<TType> join_abs(const bivector<TType>& v1, const vector<TType>& v2)
    {
    trivector<TType> res;
    res.e0e1e2 = fabs(v1.e0e1 * v2.e2) + fabs(v1.e0e2 * v2.e1) + fabs(v1.e1e2 * v2.e0);
    res.e0e1e3 = fabs(v1.e0e1 * v2.e3) + fabs(v1.e0e3 * v2.e1) + fabs(v1.e1e3 * v2.e0);
    res.e0e2e3 = fabs(v1.e0e2 * v2.e3) + fabs(v1.e0e3 * v2.e2) + fabs(v1.e2e3 * v2.e0);
    res.e1e2e3 = fabs(v1.e1e2 * v2.e3) + fabs(v1.e1e3 * v2.e2) + fabs(v1.e2e3 * v2.e1);
    return res;
    }

  template <class TType>
  trivector<TType> join_abs(const vector<TType>& v1, const bivector<TType>& v2)
    {
    return join_abs(v2, v1);
    }

  /*
  Duality of meet and join:
  (A /cap B)^* = B^* \cup A^*
  where /cap is the intersection (meet)
  /cup is the union (join)
  A^* is the dual of A
  */

  template <class TType>
  bivector<TType> meet(const trivector<TType>& v1, const trivector<TType>& v2)
    {
    auto v1_dual = dual(v1);
    auto v2_dual = dual(v2);
    auto j = join(v1_dual, v2_dual);
    return undo_dual(j);
    }

  template <class TType>
  vector<TType> meet(const bivector<TType>& v1, const trivector<TType>& v2)
    {
    auto v1_dual = dual(v1);
    auto v2_dual = dual(v2);
    auto j = join(v1_dual, v2_dual);
    return undo_dual(j);
    }

  template <class TType>
  vector<TType> meet(const trivector<TType>& v1, const bivector<TType>& v2)
    {
    auto v1_dual = dual(v1);
    auto v2_dual = dual(v2);
    auto j = join(v1_dual, v2_dual);
    return undo_dual(j);
    }

  template <class TType>
  bivector<TType> meet_abs(const trivector<TType>& v1, const trivector<TType>& v2)
    {
    auto v1_dual = dual(v1);
    auto v2_dual = dual(v2);
    auto j = join_abs(v1_dual, v2_dual);
    return abs(undo_dual(j));
    }

  template <class TType>
  vector<TType> meet_abs(const bivector<TType>& v1, const trivector<TType>& v2)
    {
    auto v1_dual = dual(v1);
    auto v2_dual = dual(v2);
    auto j = join_abs(v1_dual, v2_dual);
    return abs(undo_dual(j));
    }

  template <class TType>
  vector<TType> meet_abs(const trivector<TType>& v1, const bivector<TType>& v2)
    {
    auto v1_dual = dual(v1);
    auto v2_dual = dual(v2);
    auto j = join_abs(v1_dual, v2_dual);
    return abs(undo_dual(j));
    }

  /*
  output value: 1: true
               -1: false
                0: degenerate
  */
  template <class TType>
  int triangle_edge_intersect(const vector<TType>& v0, const vector<TType>& v1, const vector<TType>& v2, const vector<TType>& e0, const vector<TType>& e1)
    {
    trivector<TType> tria = join(join(v0, v1), v2);
    bivector<TType> edge = join(e0, e1);

    vector<TType> isct = meet(edge, tria);

    if (isct.e3 == 0)
      return 0;

    if (isct.e3 < TType(0)) // negative homogeneous coordinate
      isct = -isct;

    trivector<TType> tr0 = join(join(isct, v1), v2);
    trivector<TType> tr1 = join(join(v0, isct), v2);
    trivector<TType> tr2 = join(join(v0, v1), isct);

    bivector<TType> edge0 = join(isct, e1);
    bivector<TType> edge1 = join(e0, isct);

    auto intr0 = inner(tria, tr0);
    if (intr0 < TType(0))
      return -1;
    auto intr1 = inner(tria, tr1);
    if (intr1 < TType(0))
      return -1;
    auto intr2 = inner(tria, tr2);
    if (intr2 < TType(0))
      return -1;
    auto ined0 = inner(edge, edge0);
    if (ined0 < TType(0))
      return -1;
    auto ined1 = inner(edge, edge1);
    if (ined1 < TType(0))
      return -1;

    if (intr0 == 0 || intr1 == 0 || intr2 == 0 || ined0 == 0 || ined1 == 0)
      return 0;

    return 1;
    }

  template <class TType>
  vector<TType> triangle_edge_intersection_coordinate(const vector<TType>& v0, const vector<TType>& v1, const vector<TType>& v2, const vector<TType>& e0, const vector<TType>& e1)
    {
    trivector<TType> tria = join(join(v0, v1), v2);
    bivector<TType> edge = join(e0, e1);

    vector<TType> isct = meet(edge, tria);
    return isct;
    }

  /*
  output value: 1: true
               -1: false
                0: degenerate
  */
  template <class TType>
  int triangle_triangle_triangle_intersect(const vector<TType>& u0, const vector<TType>& u1, const vector<TType>& u2,
    const vector<TType>& v0, const vector<TType>& v1, const vector<TType>& v2,
    const vector<TType>& w0, const vector<TType>& w1, const vector<TType>& w2)
    {
    trivector<TType> tria[3];
    tria[0] = join(join(u0, u1), u2);
    tria[1] = join(join(v0, v1), v2);
    tria[2] = join(join(w0, w1), w2);

    vector<TType> isct = meet(meet(tria[0], tria[1]), tria[2]);
    if (isct.e3 == 0)
      return 0;
    if (isct.e3 < TType(0)) // negative homogeneous coordinate
      isct = -isct;

    auto intr00 = inner(tria[0], join(join(isct, u1), u2));
    if (intr00 < TType(0))
      return -1;
    auto intr01 = inner(tria[0], join(join(u0, isct), u2));
    if (intr01 < TType(0))
      return -1;
    auto intr02 = inner(tria[0], join(join(u0, u1), isct));
    if (intr02 < TType(0))
      return -1;

    auto intr10 = inner(tria[1], join(join(isct, v1), v2));
    if (intr10 < TType(0))
      return -1;
    auto intr11 = inner(tria[1], join(join(v0, isct), v2));
    if (intr11 < TType(0))
      return -1;
    auto intr12 = inner(tria[1], join(join(v0, v1), isct));
    if (intr12 < TType(0))
      return -1;

    auto intr20 = inner(tria[2], join(join(isct, w1), w2));
    if (intr20 < TType(0))
      return -1;
    auto intr21 = inner(tria[2], join(join(w0, isct), w2));
    if (intr21 < TType(0))
      return -1;
    auto intr22 = inner(tria[2], join(join(w0, w1), isct));
    if (intr22 < TType(0))
      return -1;

    if (intr00 == 0 || intr01 == 0 || intr02 == 0 || intr10 == 0 || intr11 == 0 || intr12 == 0 || intr20 == 0 || intr21 == 0 || intr22 == 0)
      return 0;

    return 1;
    }

  template <class TType>
  vector<TType> triangle_triangle_triangle_intersection_coordinate(const vector<TType>& u0, const vector<TType>& u1, const vector<TType>& u2,
    const vector<TType>& v0, const vector<TType>& v1, const vector<TType>& v2,
    const vector<TType>& w0, const vector<TType>& w1, const vector<TType>& w2)
    {
    trivector<TType> tria[3];
    tria[0] = join(join(u0, u1), u2);
    tria[1] = join(join(v0, v1), v2);
    tria[2] = join(join(w0, w1), w2);

    vector<TType> isct = meet(meet(tria[0], tria[1]), tria[2]);
    return isct;
    }

  /*
  output value: 1: true
               -1: false
                0: uncertain, more precision is needed
  */
  template <class TType>
  int triangle_edge_intersect_floating_point_filter(const vector<TType>& v0, const vector<TType>& v1, const vector<TType>& v2, const vector<TType>& e0, const vector<TType>& e1)
    {
    trivector<TType> tria = join(join(v0, v1), v2);
    bivector<TType> edge = join(e0, e1);

    trivector<TType> tria_abs = join_abs(join_abs(v0, v1), v2);
    bivector<TType> edge_abs = join_abs(e0, e1);

    vector<TType> isct = meet(edge, tria);
    vector<TType> isct_abs = meet(edge_abs, tria_abs);

    if (fabs(isct.e3) <= isct_abs.e3*11.0*std::numeric_limits<TType>::epsilon())
      return 0;

    if (isct.e3 < 0.0) // negative homogeneous coordinate
      isct = -isct;

    trivector<TType> tr0 = join(join(isct, v1), v2);
    trivector<TType> tr1 = join(join(v0, isct), v2);
    trivector<TType> tr2 = join(join(v0, v1), isct);

    bivector<TType> edge0 = join(isct, e1);
    bivector<TType> edge1 = join(e0, isct);

    trivector<TType> tr0_abs = join_abs(join_abs(isct_abs, v1), v2);
    trivector<TType> tr1_abs = join_abs(join_abs(v0, isct_abs), v2);
    trivector<TType> tr2_abs = join_abs(join_abs(v0, v1), isct_abs);

    bivector<TType> edge0_abs = join_abs(isct_abs, e1);
    bivector<TType> edge1_abs = join_abs(e0, isct_abs);

    bool uncertain = false;
    auto dot = inner(tria, tr0);
    auto dot_abs = inner(tria_abs, tr0_abs);
    bool outside = dot < 0.0;
    bool reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)   
      uncertain = true;

    dot = inner(tria, tr1);
    dot_abs = inner(tria_abs, tr1_abs);
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    dot = inner(tria, tr2);
    dot_abs = inner(tria_abs, tr2_abs);
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    dot = inner(edge, edge0);
    dot_abs = inner(edge_abs, edge0_abs);
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*21.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    dot = inner(edge, edge1);
    dot_abs = inner(edge_abs, edge1_abs);
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*21.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    if (uncertain)
      return 0;    
    return 1; // i.e. true
    }
  
  /*
  output value: 1: true
               -1: false
                0: uncertain, more precision is needed
  */
  template <class TType>
  int triangle_triangle_triangle_intersect_floating_point_filter(const vector<TType>& u0, const vector<TType>& u1, const vector<TType>& u2,
    const vector<TType>& v0, const vector<TType>& v1, const vector<TType>& v2,
    const vector<TType>& w0, const vector<TType>& w1, const vector<TType>& w2)
    {
    trivector<TType> tria[3], tria_abs[3];
    tria[0] = join(join(u0, u1), u2);
    tria[1] = join(join(v0, v1), v2);
    tria[2] = join(join(w0, w1), w2);
    tria_abs[0] = join_abs(join_abs(u0, u1), u2);
    tria_abs[1] = join_abs(join_abs(v0, v1), v2);
    tria_abs[2] = join_abs(join_abs(w0, w1), w2);

    vector<TType> isct = meet(meet(tria[0], tria[1]), tria[2]);
    vector<TType> isct_abs = meet_abs(meet_abs(tria_abs[0], tria_abs[1]), tria_abs[2]);
    if (fabs(isct.e3) <= isct_abs.e3*21.0*std::numeric_limits<TType>::epsilon())
      return 0;

    if (isct.e3 < 0.0) // negative homogeneous coordinate
      isct = -isct;

    bool uncertain = false;
    auto dot = inner(tria[0], join(join(isct, u1), u2));
    auto dot_abs = inner(tria_abs[0], join_abs(join_abs(isct_abs, u1), u2));
    bool outside = dot < 0.0;
    bool reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;
    
    dot = inner(tria[0], join(join(u0, isct), u2));
    dot_abs = inner(tria_abs[0], join_abs(join_abs(u0, isct_abs), u2));
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;
    
    dot = inner(tria[0], join(join(u0, u1), isct));
    dot_abs = inner(tria_abs[0], join_abs(join_abs(u0, u1), isct));
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    dot = inner(tria[1], join(join(isct, v1), v2));
    dot_abs = inner(tria_abs[1], join_abs(join_abs(isct_abs, v1), v2));
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    dot = inner(tria[1], join(join(v0, isct), v2));
    dot_abs = inner(tria_abs[1], join_abs(join_abs(v0, isct_abs), v2));
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    dot = inner(tria[1], join(join(v0, v1), isct));
    dot_abs = inner(tria_abs[1], join_abs(join_abs(v0, v1), isct));
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    dot = inner(tria[2], join(join(isct, w1), w2));
    dot_abs = inner(tria_abs[2], join_abs(join_abs(isct_abs, w1), w2));
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    dot = inner(tria[2], join(join(w0, isct), w2));
    dot_abs = inner(tria_abs[2], join_abs(join_abs(w0, isct_abs), w2));
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    dot = inner(tria[2], join(join(w0, w1), isct));
    dot_abs = inner(tria_abs[2], join_abs(join_abs(w0, w2), isct));
    outside = dot < 0.0;
    reliable = fabs(dot) > dot_abs*25.0*std::numeric_limits<TType>::epsilon();
    if (reliable && outside)
      return -1; // i.e. false
    if (!reliable)
      uncertain = true;

    if (uncertain)
      return 0;
    return 1; // i.e. true
    }
  }