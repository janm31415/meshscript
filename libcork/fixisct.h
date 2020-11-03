#pragma once

#include "fixint.h"
#include <stdint.h>
#include "ga4.h"

template <int N>
inline void dual(ga::vector<typename fixint::bit_int<N>::type>& out, const ga::trivector<typename fixint::bit_int<N>::type>& v)
  {
  out.e0 = v.e1e2e3;
  out.e1 = neg(v.e0e2e3);
  out.e2 = v.e0e1e3;
  out.e3 = neg(v.e0e1e2);
  }

template <int N>
inline void dual(ga::bivector<typename fixint::bit_int<N>::type>& out, const ga::bivector<typename fixint::bit_int<N>::type>& v)
  {
  out.e0e1 = v.e2e3;
  out.e0e2 = neg(v.e1e3);
  out.e0e3 = v.e1e2;
  out.e1e2 = v.e0e3;
  out.e1e3 = neg(v.e0e2);
  out.e2e3 = v.e0e1;
  }

template <int N>
inline void dual(ga::trivector<typename fixint::bit_int<N>::type>& out, const ga::vector<typename fixint::bit_int<N>::type>& v)
  {
  out.e0e1e2 = v.e3;
  out.e0e1e3 = neg(v.e2);
  out.e0e2e3 = v.e1;
  out.e1e2e3 = neg(v.e0);
  }

template <int N>
inline void undo_dual(ga::vector<typename fixint::bit_int<N>::type>& out, const ga::trivector<typename fixint::bit_int<N>::type>& v)
  {
  out.e0 = neg(v.e1e2e3);
  out.e1 = v.e0e2e3;
  out.e2 = neg(v.e0e1e3);
  out.e3 = v.e0e1e2;
  }

template <int N>
inline void undo_dual(ga::bivector<typename fixint::bit_int<N>::type>& out, const ga::bivector<typename fixint::bit_int<N>::type>& v)
  {
  out.e0e1 = v.e2e3;
  out.e0e2 = neg(v.e1e3);
  out.e0e3 = v.e1e2;
  out.e1e2 = v.e0e3;
  out.e1e3 = neg(v.e0e2);
  out.e2e3 = v.e0e1;
  }

template <int N>
inline void undo_dual(ga::trivector<typename fixint::bit_int<N>::type>& out, const ga::vector<typename fixint::bit_int<N>::type>& v)
  {
  out.e0e1e2 = neg(v.e3);
  out.e0e1e3 = v.e2;
  out.e0e2e3 = neg(v.e1);
  out.e1e2e3 = v.e0;
  }

template <int Nlhs, int Nrhs>
inline void join(ga::bivector<typename fixint::bit_int<Nlhs + Nrhs + 1>::type>& out, const ga::vector<typename fixint::bit_int<Nlhs>::type>& v1, const ga::vector<typename fixint::bit_int<Nrhs>::type>& v2)
  {
  /*
  out.e0e1 = (v1.e0 * v2.e1) - (v2.e0 * v1.e1);
  out.e0e2 = (v1.e0 * v2.e2) - (v2.e0 * v1.e2);
  out.e0e3 = (v1.e0 * v2.e3) - (v2.e0 * v1.e3);
  out.e1e2 = (v1.e1 * v2.e2) - (v2.e1 * v1.e2);
  out.e1e3 = (v1.e1 * v2.e3) - (v2.e1 * v1.e3);
  out.e2e3 = (v1.e2 * v2.e3) - (v2.e2 * v1.e3);
  */
  using namespace fixint;
  typename bit_int<Nlhs + Nrhs>::type a;
  typename bit_int<Nlhs + Nrhs>::type b;
  mul(a, v1.e0, v2.e1);
  mul(b, v2.e0, v1.e1);
  sub(out.e0e1, a, b);

  mul(a, v1.e0, v2.e2);
  mul(b, v2.e0, v1.e2);
  sub(out.e0e2, a, b);

  mul(a, v1.e0, v2.e3);
  mul(b, v2.e0, v1.e3);
  sub(out.e0e3, a, b);

  mul(a, v1.e1, v2.e2);
  mul(b, v2.e1, v1.e2);
  sub(out.e1e2, a, b);

  mul(a, v1.e1, v2.e3);
  mul(b, v2.e1, v1.e3);
  sub(out.e1e3, a, b);

  mul(a, v1.e2, v2.e3);
  mul(b, v2.e2, v1.e3);
  sub(out.e2e3, a, b);
  }

template <int Nlhs, int Nrhs>
inline void join(ga::trivector<typename fixint::bit_int<Nlhs + Nrhs + 2>::type>& out, const ga::bivector<typename fixint::bit_int<Nlhs>::type>& v1, const ga::vector<typename fixint::bit_int<Nrhs>::type>& v2)
  {
  /*
  out.e0e1e2 = (v1.e0e1 * v2.e2) - (v1.e0e2 * v2.e1) + (v1.e1e2 * v2.e0);
  out.e0e1e3 = (v1.e0e1 * v2.e3) - (v1.e0e3 * v2.e1) + (v1.e1e3 * v2.e0);
  out.e0e2e3 = (v1.e0e2 * v2.e3) - (v1.e0e3 * v2.e2) + (v1.e2e3 * v2.e0);
  out.e1e2e3 = (v1.e1e2 * v2.e3) - (v1.e1e3 * v2.e2) + (v1.e2e3 * v2.e1);
  */
  using namespace fixint;
  typename bit_int<Nlhs + Nrhs>::type a;
  typename bit_int<Nlhs + Nrhs>::type b;
  typename bit_int<Nlhs + Nrhs>::type c;
  typename bit_int<Nlhs + Nrhs + 1>::type x;
  mul(a, v1.e0e1, v2.e2);
  mul(b, v1.e0e2, v2.e1);
  mul(c, v1.e1e2, v2.e0);
  sub(x, a, b);
  add(out.e0e1e2, x, c);

  mul(a, v1.e0e1, v2.e3);
  mul(b, v1.e0e3, v2.e1);
  mul(c, v1.e1e3, v2.e0);
  sub(x, a, b);
  add(out.e0e1e3, x, c);

  mul(a, v1.e0e2, v2.e3);
  mul(b, v1.e0e3, v2.e2);
  mul(c, v1.e2e3, v2.e0);
  sub(x, a, b);
  add(out.e0e2e3, x, c);

  mul(a, v1.e1e2, v2.e3);
  mul(b, v1.e1e3, v2.e2);
  mul(c, v1.e2e3, v2.e1);
  sub(x, a, b);
  add(out.e1e2e3, x, c);
  }

template <int Nlhs, int Nrhs>
inline void meet(ga::vector<typename fixint::bit_int<Nlhs + Nrhs + 2>::type>& out, const ga::bivector<typename fixint::bit_int<Nlhs>::type>& v1, const ga::trivector<typename fixint::bit_int<Nrhs>::type>& v2)
  {
  /*
  auto v1_dual = dual(v1);
  auto v2_dual = dual(v2);
  auto j = join(v1_dual, v2_dual);
  return undo_dual(j);
  */
  using namespace fixint;
  ga::trivector<typename bit_int<Nlhs + Nrhs + 2>::type> out_dual;
  ga::bivector<typename bit_int<Nlhs>::type> v1_dual;
  ga::vector<typename bit_int<Nrhs>::type> v2_dual;
  dual<Nlhs>(v1_dual, v1);
  dual<Nrhs>(v2_dual, v2);
  join<Nlhs, Nrhs>(out_dual, v1_dual, v2_dual);
  undo_dual<Nlhs + Nrhs + 2>(out, out_dual);
  }

template <int Nlhs, int Nrhs>
inline void meet(ga::bivector<typename fixint::bit_int<Nlhs + Nrhs + 1>::type>& out, const ga::trivector<typename fixint::bit_int<Nlhs>::type>& v1, const ga::trivector<typename fixint::bit_int<Nrhs>::type>& v2)
  {
  /*
    auto v1_dual = dual(v1);
    auto v2_dual = dual(v2);
    auto j = join(v1_dual, v2_dual);
    return undo_dual(j);
  */
  using namespace fixint;
  ga::bivector<typename bit_int<Nlhs + Nrhs + 1>::type> out_dual;
  ga::vector<typename bit_int<Nlhs>::type> v1_dual;
  ga::vector<typename bit_int<Nrhs>::type> v2_dual;
  dual<Nlhs>(v1_dual, v1);
  dual<Nrhs>(v2_dual, v2);
  join<Nlhs, Nrhs>(out_dual, v1_dual, v2_dual);
  undo_dual<Nlhs + Nrhs + 1>(out, out_dual);
  }

template <int Nlhs, int Nrhs>
inline void inner(typename fixint::bit_int<Nlhs + Nrhs + 2>::type& out, const ga::trivector<typename fixint::bit_int<Nlhs>::type>& v1, const ga::trivector<typename fixint::bit_int<Nrhs>::type>& v2)
  {
  //return v1.e0e1e2*v2.e0e1e2 + v1.e0e1e3*v2.e0e1e3 + v1.e0e2e3*v2.e0e2e3 + v1.e1e2e3*v2.e1e2e3;
  using namespace fixint;
  typename bit_int<Nlhs + Nrhs>::type m0, m1, m2, m3;
  mul(m0, v1.e0e1e2, v2.e0e1e2);
  mul(m1, v1.e0e1e3, v2.e0e1e3);
  mul(m2, v1.e0e2e3, v2.e0e2e3);
  mul(m3, v1.e1e2e3, v2.e1e2e3);
  typename bit_int<Nlhs + Nrhs + 1>::type a, b;
  add(a, m0, m1);
  add(b, m2, m3);
  add(out, a, b);
  }

template <int Nlhs, int Nrhs>
inline void inner(typename fixint::bit_int<Nlhs + Nrhs + 3>::type& out, const ga::bivector<typename fixint::bit_int<Nlhs>::type>& v1, const ga::bivector<typename fixint::bit_int<Nrhs>::type>& v2)
  {
  // return v1.e0e1*v2.e0e1 + v1.e0e2*v2.e0e2 + v1.e0e3*v2.e0e3 + v1.e1e2*v2.e1e2 + v1.e1e3*v2.e1e3 + v1.e2e3*v2.e2e3;
  using namespace fixint;
  typename bit_int<Nlhs + Nrhs>::type m0, m1, m2, m3, m4, m5;
  mul(m0, v1.e0e1, v2.e0e1);
  mul(m1, v1.e0e2, v2.e0e2);
  mul(m2, v1.e0e3, v2.e0e3);
  mul(m3, v1.e1e2, v2.e1e2);
  mul(m4, v1.e1e3, v2.e1e3);
  mul(m5, v1.e2e3, v2.e2e3);
  typename bit_int<Nlhs + Nrhs + 1>::type a, b, c;
  add(a, m0, m1);
  add(b, m2, m3);
  add(c, m4, m5);
  typename bit_int<Nlhs + Nrhs + 2>::type d;
  add(d, a, b);
  add(out, d, c);
  }

inline int triangle_edge_intersect_fixed(const ga::vector<fixint::limb_int<1>>& v0, const ga::vector<fixint::limb_int<1>>& v1, const ga::vector<fixint::limb_int<1>>& v2, const ga::vector<fixint::limb_int<1>>& e0, const ga::vector<fixint::limb_int<1>>& e1)
  {
  const static int IN_BITS = 31;
  const static int LINE_BITS = 2 * IN_BITS + 1;
  const static int TRI_BITS = LINE_BITS + IN_BITS + 2;
  const static int ISCT_BITS = TRI_BITS + LINE_BITS + 2;
  const static int LINE_A_BITS = ISCT_BITS + IN_BITS + 1;
  const static int TRI_A_BITS = LINE_A_BITS + IN_BITS + 2;
  const static int INNER_LINE_BITS = LINE_A_BITS + LINE_BITS + 3;
  const static int INNER_TRI_BITS = TRI_A_BITS + TRI_BITS + 2;

  ga::bivector<typename fixint::bit_int<LINE_BITS>::type> temp;
  join<IN_BITS, IN_BITS>(temp, v0, v1);
  ga::trivector<typename fixint::bit_int<TRI_BITS>::type> tria;
  join<LINE_BITS, IN_BITS>(tria, temp, v2);
  ga::bivector<typename fixint::bit_int<LINE_BITS>::type> edge;
  join<IN_BITS, IN_BITS>(edge, e0, e1);
  ga::vector<typename fixint::bit_int<ISCT_BITS>::type> isct;
  meet<LINE_BITS, TRI_BITS>(isct, edge, tria);

  if (fixint::is_zero(isct.e3))
    return 0;

  if (fixint::is_negative(isct.e3))
    {
    isct.e0 = fixint::neg(isct.e0);
    isct.e1 = fixint::neg(isct.e1);
    isct.e2 = fixint::neg(isct.e2);
    isct.e3 = fixint::neg(isct.e3);
    }

  ga::bivector<typename fixint::bit_int<LINE_A_BITS>::type> temp2;
  ga::trivector<typename fixint::bit_int<TRI_A_BITS>::type> tr0, tr1, tr2;
  join<ISCT_BITS, IN_BITS>(temp2, isct, v1);
  join<LINE_A_BITS, IN_BITS>(tr0, temp2, v2);
  join<IN_BITS, ISCT_BITS>(temp2, v0, isct);
  join<LINE_A_BITS, IN_BITS>(tr1, temp2, v2);
  join<LINE_BITS, ISCT_BITS>(tr2, temp, isct);

  ga::bivector<typename fixint::bit_int<LINE_A_BITS>::type> edge0, edge1;
  join<ISCT_BITS, IN_BITS>(edge0, isct, e1);
  join<IN_BITS, ISCT_BITS>(edge1, e0, isct);

  typename fixint::bit_int<INNER_TRI_BITS>::type tri_dot0, tri_dot1, tri_dot2;
  inner<TRI_BITS, TRI_A_BITS>(tri_dot0, tria, tr0);
  if (fixint::is_negative(tri_dot0))
    return -1;
  inner<TRI_BITS, TRI_A_BITS>(tri_dot1, tria, tr1);
  if (fixint::is_negative(tri_dot1))
    return -1;
  inner<TRI_BITS, TRI_A_BITS>(tri_dot2, tria, tr2);
  if (fixint::is_negative(tri_dot2))
    return -1;
  typename fixint::bit_int<INNER_LINE_BITS>::type line_dot0, line_dot1;
  inner<LINE_BITS, LINE_A_BITS>(line_dot0, edge, edge0);
  if (fixint::is_negative(line_dot0))
    return -1;
  inner<LINE_BITS, LINE_A_BITS>(line_dot1, edge, edge1);
  if (fixint::is_negative(line_dot1))
    return -1;
  if (fixint::is_zero(tri_dot0) || fixint::is_zero(tri_dot1) || fixint::is_zero(tri_dot2) || fixint::is_zero(line_dot0) || fixint::is_zero(line_dot1))
    return 0;
  return 1;
  }


inline int triangle_edge_intersect_fixed(const ga::vector<fixint::limb_int<2>>& v0, const ga::vector<fixint::limb_int<2>>& v1, const ga::vector<fixint::limb_int<2>>& v2, const ga::vector<fixint::limb_int<2>>& e0, const ga::vector<fixint::limb_int<2>>& e1)
  {
  const static int IN_BITS = 63;
  const static int LINE_BITS = 2 * IN_BITS + 1;
  const static int TRI_BITS = LINE_BITS + IN_BITS + 2;
  const static int ISCT_BITS = TRI_BITS + LINE_BITS + 2;
  const static int LINE_A_BITS = ISCT_BITS + IN_BITS + 1;
  const static int TRI_A_BITS = LINE_A_BITS + IN_BITS + 2;
  const static int INNER_LINE_BITS = LINE_A_BITS + LINE_BITS + 3;
  const static int INNER_TRI_BITS = TRI_A_BITS + TRI_BITS + 2;

  ga::bivector<typename fixint::bit_int<LINE_BITS>::type> temp;
  join<IN_BITS, IN_BITS>(temp, v0, v1);
  ga::trivector<typename fixint::bit_int<TRI_BITS>::type> tria;
  join<LINE_BITS, IN_BITS>(tria, temp, v2);
  ga::bivector<typename fixint::bit_int<LINE_BITS>::type> edge;
  join<IN_BITS, IN_BITS>(edge, e0, e1);
  ga::vector<typename fixint::bit_int<ISCT_BITS>::type> isct;
  meet<LINE_BITS, TRI_BITS>(isct, edge, tria);

  if (fixint::is_zero(isct.e3))
    return 0;

  if (fixint::is_negative(isct.e3))
    {
    isct.e0 = fixint::neg(isct.e0);
    isct.e1 = fixint::neg(isct.e1);
    isct.e2 = fixint::neg(isct.e2);
    isct.e3 = fixint::neg(isct.e3);
    }

  ga::bivector<typename fixint::bit_int<LINE_A_BITS>::type> temp2;
  ga::trivector<typename fixint::bit_int<TRI_A_BITS>::type> tr0, tr1, tr2;
  join<ISCT_BITS, IN_BITS>(temp2, isct, v1);
  join<LINE_A_BITS, IN_BITS>(tr0, temp2, v2);
  join<IN_BITS, ISCT_BITS>(temp2, v0, isct);
  join<LINE_A_BITS, IN_BITS>(tr1, temp2, v2);
  join<LINE_BITS, ISCT_BITS>(tr2, temp, isct);

  ga::bivector<typename fixint::bit_int<LINE_A_BITS>::type> edge0, edge1;
  join<ISCT_BITS, IN_BITS>(edge0, isct, e1);
  join<IN_BITS, ISCT_BITS>(edge1, e0, isct);

  typename fixint::bit_int<INNER_TRI_BITS>::type tri_dot0, tri_dot1, tri_dot2;
  inner<TRI_BITS, TRI_A_BITS>(tri_dot0, tria, tr0);
  if (fixint::is_negative(tri_dot0))
    return -1;
  inner<TRI_BITS, TRI_A_BITS>(tri_dot1, tria, tr1);
  if (fixint::is_negative(tri_dot1))
    return -1;
  inner<TRI_BITS, TRI_A_BITS>(tri_dot2, tria, tr2);
  if (fixint::is_negative(tri_dot2))
    return -1;
  typename fixint::bit_int<INNER_LINE_BITS>::type line_dot0, line_dot1;
  inner<LINE_BITS, LINE_A_BITS>(line_dot0, edge, edge0);
  if (fixint::is_negative(line_dot0))
    return -1;
  inner<LINE_BITS, LINE_A_BITS>(line_dot1, edge, edge1);
  if (fixint::is_negative(line_dot1))
    return -1;
  if (fixint::is_zero(tri_dot0) || fixint::is_zero(tri_dot1) || fixint::is_zero(tri_dot2) || fixint::is_zero(line_dot0) || fixint::is_zero(line_dot1))
    return 0;
  return 1;
  }


inline int triangle_edge_intersect_fixed(const ga::vector<int32_t>& v0, const ga::vector<int32_t>& v1, const ga::vector<int32_t>& v2, const ga::vector<int32_t>& e0, const ga::vector<int32_t>& e1)
  {
  ga::vector<fixint::limb_int<1>> V0, V1, V2, E0, E1;
  V0.e0 = fixint::limb_int<1>(v0.e0);
  V0.e1 = fixint::limb_int<1>(v0.e1);
  V0.e2 = fixint::limb_int<1>(v0.e2);
  V0.e3 = fixint::limb_int<1>(v0.e3);
  V1.e0 = fixint::limb_int<1>(v1.e0);
  V1.e1 = fixint::limb_int<1>(v1.e1);
  V1.e2 = fixint::limb_int<1>(v1.e2);
  V1.e3 = fixint::limb_int<1>(v1.e3);
  V2.e0 = fixint::limb_int<1>(v2.e0);
  V2.e1 = fixint::limb_int<1>(v2.e1);
  V2.e2 = fixint::limb_int<1>(v2.e2);
  V2.e3 = fixint::limb_int<1>(v2.e3);
  E0.e0 = fixint::limb_int<1>(e0.e0);
  E0.e1 = fixint::limb_int<1>(e0.e1);
  E0.e2 = fixint::limb_int<1>(e0.e2);
  E0.e3 = fixint::limb_int<1>(e0.e3);
  E1.e0 = fixint::limb_int<1>(e1.e0);
  E1.e1 = fixint::limb_int<1>(e1.e1);
  E1.e2 = fixint::limb_int<1>(e1.e2);
  E1.e3 = fixint::limb_int<1>(e1.e3);
  return triangle_edge_intersect_fixed(V0, V1, V2, E0, E1);
  }

inline int triangle_triangle_triangle_intersect_fixed(const ga::vector<fixint::limb_int<1>>& u0, const ga::vector<fixint::limb_int<1>>& u1, const ga::vector<fixint::limb_int<1>>& u2, const ga::vector<fixint::limb_int<1>>& v0, const ga::vector<fixint::limb_int<1>>& v1, const ga::vector<fixint::limb_int<1>>& v2, const ga::vector<fixint::limb_int<1>>& w0, const ga::vector<fixint::limb_int<1>>& w1, const ga::vector<fixint::limb_int<1>>& w2)
  {
  const static int IN_BITS = 31;
  const static int EXT2_UP_BITS = 2 * IN_BITS + 1;
  const static int EXT3_UP_BITS = EXT2_UP_BITS + IN_BITS + 2;
  const static int EXT2_DN_BITS = 2 * EXT3_UP_BITS + 1;
  const static int ISCT_BITS = EXT2_DN_BITS + EXT3_UP_BITS + 2;
  const static int EXT2_TA_BITS = ISCT_BITS + IN_BITS + 1;
  const static int EXT3_TA_BITS = EXT2_TA_BITS + IN_BITS + 2;
  const static int INNER_BITS = EXT3_TA_BITS + EXT3_UP_BITS + 2;

  ga::bivector<typename fixint::bit_int<EXT2_UP_BITS>::type> temp;
  ga::trivector<typename fixint::bit_int<EXT3_UP_BITS>::type> tria[3];
  join<IN_BITS, IN_BITS>(temp, u0, u1);
  join<EXT2_UP_BITS, IN_BITS>(tria[0], temp, u2);
  join<IN_BITS, IN_BITS>(temp, v0, v1);
  join<EXT2_UP_BITS, IN_BITS>(tria[1], temp, v2);
  join<IN_BITS, IN_BITS>(temp, w0, w1);
  join<EXT2_UP_BITS, IN_BITS>(tria[2], temp, w2);

  ga::bivector<typename fixint::bit_int<EXT2_DN_BITS>::type> isct_edge;
  ga::vector<typename fixint::bit_int<ISCT_BITS>::type> isct;
  meet<EXT3_UP_BITS, EXT3_UP_BITS>(isct_edge, tria[0], tria[1]);
  meet<EXT2_DN_BITS, EXT3_UP_BITS>(isct, isct_edge, tria[2]);

  if (fixint::is_zero(isct.e3))
    return 0;

  if (fixint::is_negative(isct.e3))
    {
    isct.e0 = fixint::neg(isct.e0);
    isct.e1 = fixint::neg(isct.e1);
    isct.e2 = fixint::neg(isct.e2);
    isct.e3 = fixint::neg(isct.e3);
    }

  ga::trivector<typename fixint::bit_int<EXT3_TA_BITS>::type> a;
  ga::bivector<typename fixint::bit_int<EXT2_TA_BITS>::type> t;
  ga::bivector<typename fixint::bit_int<EXT2_UP_BITS>::type> t2;

  join<ISCT_BITS, IN_BITS>(t, isct, u1);
  join<EXT2_TA_BITS, IN_BITS>(a, t, u2);
  typename fixint::bit_int<INNER_BITS>::type intr00;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr00, tria[0], a);
  if (fixint::is_negative(intr00))
    return -1;

  join<IN_BITS, ISCT_BITS>(t, u0, isct);
  join<EXT2_TA_BITS, IN_BITS>(a, t, u2);
  typename fixint::bit_int<INNER_BITS>::type intr01;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr01, tria[0], a);
  if (fixint::is_negative(intr01))
    return -1;

  join<IN_BITS, IN_BITS>(t2, u0, u1);
  join<EXT2_UP_BITS, ISCT_BITS>(a, t2, isct);
  typename fixint::bit_int<INNER_BITS>::type intr02;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr02, tria[0], a);
  if (fixint::is_negative(intr02))
    return -1;

  join<ISCT_BITS, IN_BITS>(t, isct, v1);
  join<EXT2_TA_BITS, IN_BITS>(a, t, v2);
  typename fixint::bit_int<INNER_BITS>::type intr10;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr10, tria[1], a);
  if (fixint::is_negative(intr10))
    return -1;

  join<IN_BITS, ISCT_BITS>(t, v0, isct);
  join<EXT2_TA_BITS, IN_BITS>(a, t, v2);
  typename fixint::bit_int<INNER_BITS>::type intr11;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr11, tria[1], a);
  if (fixint::is_negative(intr11))
    return -1;

  join<IN_BITS, IN_BITS>(t2, v0, v1);
  join<EXT2_UP_BITS, ISCT_BITS>(a, t2, isct);
  typename fixint::bit_int<INNER_BITS>::type intr12;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr12, tria[1], a);
  if (fixint::is_negative(intr12))
    return -1;

  join<ISCT_BITS, IN_BITS>(t, isct, w1);
  join<EXT2_TA_BITS, IN_BITS>(a, t, w2);
  typename fixint::bit_int<INNER_BITS>::type intr20;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr20, tria[2], a);
  if (fixint::is_negative(intr20))
    return -1;

  join<IN_BITS, ISCT_BITS>(t, w0, isct);
  join<EXT2_TA_BITS, IN_BITS>(a, t, w2);
  typename fixint::bit_int<INNER_BITS>::type intr21;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr21, tria[2], a);
  if (fixint::is_negative(intr21))
    return -1;

  join<IN_BITS, IN_BITS>(t2, w0, w1);
  join<EXT2_UP_BITS, ISCT_BITS>(a, t2, isct);
  typename fixint::bit_int<INNER_BITS>::type intr22;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr22, tria[2], a);
  if (fixint::is_negative(intr22))
    return -1;

  if (fixint::is_zero(intr00) || fixint::is_zero(intr01) || fixint::is_zero(intr02) || fixint::is_zero(intr10) || fixint::is_zero(intr11) || fixint::is_zero(intr12) || fixint::is_zero(intr20) || fixint::is_zero(intr21) || fixint::is_zero(intr22))
    return 0;

  return 1;
  }


inline int triangle_triangle_triangle_intersect_fixed(const ga::vector<fixint::limb_int<2>>& u0, const ga::vector<fixint::limb_int<2>>& u1, const ga::vector<fixint::limb_int<2>>& u2, const ga::vector<fixint::limb_int<2>>& v0, const ga::vector<fixint::limb_int<2>>& v1, const ga::vector<fixint::limb_int<2>>& v2, const ga::vector<fixint::limb_int<2>>& w0, const ga::vector<fixint::limb_int<2>>& w1, const ga::vector<fixint::limb_int<2>>& w2)
  {
  const static int IN_BITS = 63;
  const static int EXT2_UP_BITS = 2 * IN_BITS + 1;
  const static int EXT3_UP_BITS = EXT2_UP_BITS + IN_BITS + 2;
  const static int EXT2_DN_BITS = 2 * EXT3_UP_BITS + 1;
  const static int ISCT_BITS = EXT2_DN_BITS + EXT3_UP_BITS + 2;
  const static int EXT2_TA_BITS = ISCT_BITS + IN_BITS + 1;
  const static int EXT3_TA_BITS = EXT2_TA_BITS + IN_BITS + 2;
  const static int INNER_BITS = EXT3_TA_BITS + EXT3_UP_BITS + 2;

  ga::bivector<typename fixint::bit_int<EXT2_UP_BITS>::type> temp;
  ga::trivector<typename fixint::bit_int<EXT3_UP_BITS>::type> tria[3];
  join<IN_BITS, IN_BITS>(temp, u0, u1);
  join<EXT2_UP_BITS, IN_BITS>(tria[0], temp, u2);
  join<IN_BITS, IN_BITS>(temp, v0, v1);
  join<EXT2_UP_BITS, IN_BITS>(tria[1], temp, v2);
  join<IN_BITS, IN_BITS>(temp, w0, w1);
  join<EXT2_UP_BITS, IN_BITS>(tria[2], temp, w2);

  ga::bivector<typename fixint::bit_int<EXT2_DN_BITS>::type> isct_edge;
  ga::vector<typename fixint::bit_int<ISCT_BITS>::type> isct;
  meet<EXT3_UP_BITS, EXT3_UP_BITS>(isct_edge, tria[0], tria[1]);
  meet<EXT2_DN_BITS, EXT3_UP_BITS>(isct, isct_edge, tria[2]);

  if (fixint::is_zero(isct.e3))
    return 0;

  if (fixint::is_negative(isct.e3))
    {
    isct.e0 = fixint::neg(isct.e0);
    isct.e1 = fixint::neg(isct.e1);
    isct.e2 = fixint::neg(isct.e2);
    isct.e3 = fixint::neg(isct.e3);
    }

  ga::trivector<typename fixint::bit_int<EXT3_TA_BITS>::type> a;
  ga::bivector<typename fixint::bit_int<EXT2_TA_BITS>::type> t;
  ga::bivector<typename fixint::bit_int<EXT2_UP_BITS>::type> t2;

  join<ISCT_BITS, IN_BITS>(t, isct, u1);
  join<EXT2_TA_BITS, IN_BITS>(a, t, u2);
  typename fixint::bit_int<INNER_BITS>::type intr00;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr00, tria[0], a);
  if (fixint::is_negative(intr00))
    return -1;

  join<IN_BITS, ISCT_BITS>(t, u0, isct);
  join<EXT2_TA_BITS, IN_BITS>(a, t, u2);
  typename fixint::bit_int<INNER_BITS>::type intr01;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr01, tria[0], a);
  if (fixint::is_negative(intr01))
    return -1;

  join<IN_BITS, IN_BITS>(t2, u0, u1);
  join<EXT2_UP_BITS, ISCT_BITS>(a, t2, isct);
  typename fixint::bit_int<INNER_BITS>::type intr02;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr02, tria[0], a);
  if (fixint::is_negative(intr02))
    return -1;

  join<ISCT_BITS, IN_BITS>(t, isct, v1);
  join<EXT2_TA_BITS, IN_BITS>(a, t, v2);
  typename fixint::bit_int<INNER_BITS>::type intr10;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr10, tria[1], a);
  if (fixint::is_negative(intr10))
    return -1;

  join<IN_BITS, ISCT_BITS>(t, v0, isct);
  join<EXT2_TA_BITS, IN_BITS>(a, t, v2);
  typename fixint::bit_int<INNER_BITS>::type intr11;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr11, tria[1], a);
  if (fixint::is_negative(intr11))
    return -1;

  join<IN_BITS, IN_BITS>(t2, v0, v1);
  join<EXT2_UP_BITS, ISCT_BITS>(a, t2, isct);
  typename fixint::bit_int<INNER_BITS>::type intr12;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr12, tria[1], a);
  if (fixint::is_negative(intr12))
    return -1;

  join<ISCT_BITS, IN_BITS>(t, isct, w1);
  join<EXT2_TA_BITS, IN_BITS>(a, t, w2);
  typename fixint::bit_int<INNER_BITS>::type intr20;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr20, tria[2], a);
  if (fixint::is_negative(intr20))
    return -1;

  join<IN_BITS, ISCT_BITS>(t, w0, isct);
  join<EXT2_TA_BITS, IN_BITS>(a, t, w2);
  typename fixint::bit_int<INNER_BITS>::type intr21;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr21, tria[2], a);
  if (fixint::is_negative(intr21))
    return -1;

  join<IN_BITS, IN_BITS>(t2, w0, w1);
  join<EXT2_UP_BITS, ISCT_BITS>(a, t2, isct);
  typename fixint::bit_int<INNER_BITS>::type intr22;
  inner<EXT3_UP_BITS, EXT3_TA_BITS>(intr22, tria[2], a);
  if (fixint::is_negative(intr22))
    return -1;

  if (fixint::is_zero(intr00) || fixint::is_zero(intr01) || fixint::is_zero(intr02) || fixint::is_zero(intr10) || fixint::is_zero(intr11) || fixint::is_zero(intr12) || fixint::is_zero(intr20) || fixint::is_zero(intr21) || fixint::is_zero(intr22))
    return 0;

  return 1;
  }


inline int triangle_triangle_triangle_intersect_fixed(const ga::vector<int32_t>& u0, const ga::vector<int32_t>& u1, const ga::vector<int32_t>& u2, const ga::vector<int32_t>& v0, const ga::vector<int32_t>& v1, const ga::vector<int32_t>& v2, const ga::vector<int32_t>& w0, const ga::vector<int32_t>& w1, const ga::vector<int32_t>& w2)
  {
  ga::vector<fixint::limb_int<1>> U0, U1, U2, V0, V1, V2, W0, W1, W2;
  U0.e0 = fixint::limb_int<1>(u0.e0);
  U0.e1 = fixint::limb_int<1>(u0.e1);
  U0.e2 = fixint::limb_int<1>(u0.e2);
  U0.e3 = fixint::limb_int<1>(u0.e3);
  U1.e0 = fixint::limb_int<1>(u1.e0);
  U1.e1 = fixint::limb_int<1>(u1.e1);
  U1.e2 = fixint::limb_int<1>(u1.e2);
  U1.e3 = fixint::limb_int<1>(u1.e3);
  U2.e0 = fixint::limb_int<1>(u2.e0);
  U2.e1 = fixint::limb_int<1>(u2.e1);
  U2.e2 = fixint::limb_int<1>(u2.e2);
  U2.e3 = fixint::limb_int<1>(u2.e3);
  V0.e0 = fixint::limb_int<1>(v0.e0);
  V0.e1 = fixint::limb_int<1>(v0.e1);
  V0.e2 = fixint::limb_int<1>(v0.e2);
  V0.e3 = fixint::limb_int<1>(v0.e3);
  V1.e0 = fixint::limb_int<1>(v1.e0);
  V1.e1 = fixint::limb_int<1>(v1.e1);
  V1.e2 = fixint::limb_int<1>(v1.e2);
  V1.e3 = fixint::limb_int<1>(v1.e3);
  V2.e0 = fixint::limb_int<1>(v2.e0);
  V2.e1 = fixint::limb_int<1>(v2.e1);
  V2.e2 = fixint::limb_int<1>(v2.e2);
  V2.e3 = fixint::limb_int<1>(v2.e3);
  W0.e0 = fixint::limb_int<1>(w0.e0);
  W0.e1 = fixint::limb_int<1>(w0.e1);
  W0.e2 = fixint::limb_int<1>(w0.e2);
  W0.e3 = fixint::limb_int<1>(w0.e3);
  W1.e0 = fixint::limb_int<1>(w1.e0);
  W1.e1 = fixint::limb_int<1>(w1.e1);
  W1.e2 = fixint::limb_int<1>(w1.e2);
  W1.e3 = fixint::limb_int<1>(w1.e3);
  W2.e0 = fixint::limb_int<1>(w2.e0);
  W2.e1 = fixint::limb_int<1>(w2.e1);
  W2.e2 = fixint::limb_int<1>(w2.e2);
  W2.e3 = fixint::limb_int<1>(w2.e3);
  return triangle_triangle_triangle_intersect_fixed(U0, U1, U2, V0, V1, V2, W0, W1, W2);
  }

inline ga::vector<typename fixint::bit_int<161>::type> triangle_edge_intersection_coordinate_fixed(const ga::vector<fixint::limb_int<1>>& v0, const ga::vector<fixint::limb_int<1>>& v1, const ga::vector<fixint::limb_int<1>>& v2, const ga::vector<fixint::limb_int<1>>& e0, const ga::vector<fixint::limb_int<1>>& e1)
  {
  const static int IN_BITS = 31;
  const static int LINE_BITS = 2 * IN_BITS + 1; //63
  const static int TRI_BITS = LINE_BITS + IN_BITS + 2; //96
  const static int ISCT_BITS = TRI_BITS + LINE_BITS + 2; //161  

  ga::bivector<typename fixint::bit_int<LINE_BITS>::type> temp;
  join<IN_BITS, IN_BITS>(temp, v0, v1);
  ga::trivector<typename fixint::bit_int<TRI_BITS>::type> tria;
  join<LINE_BITS, IN_BITS>(tria, temp, v2);
  ga::bivector<typename fixint::bit_int<LINE_BITS>::type> edge;
  join<IN_BITS, IN_BITS>(edge, e0, e1);
  ga::vector<typename fixint::bit_int<ISCT_BITS>::type> isct;
  meet<LINE_BITS, TRI_BITS>(isct, edge, tria);
  return isct;
  }

inline ga::vector<typename fixint::bit_int<321>::type> triangle_edge_intersection_coordinate_fixed(const ga::vector<fixint::limb_int<2>>& v0, const ga::vector<fixint::limb_int<2>>& v1, const ga::vector<fixint::limb_int<2>>& v2, const ga::vector<fixint::limb_int<2>>& e0, const ga::vector<fixint::limb_int<2>>& e1)
  {
  const static int IN_BITS = 63;
  const static int LINE_BITS = 2 * IN_BITS + 1; //127
  const static int TRI_BITS = LINE_BITS + IN_BITS + 2; //192
  const static int ISCT_BITS = TRI_BITS + LINE_BITS + 2; //321  

  ga::bivector<typename fixint::bit_int<LINE_BITS>::type> temp;
  join<IN_BITS, IN_BITS>(temp, v0, v1);
  ga::trivector<typename fixint::bit_int<TRI_BITS>::type> tria;
  join<LINE_BITS, IN_BITS>(tria, temp, v2);
  ga::bivector<typename fixint::bit_int<LINE_BITS>::type> edge;
  join<IN_BITS, IN_BITS>(edge, e0, e1);
  ga::vector<typename fixint::bit_int<ISCT_BITS>::type> isct;
  meet<LINE_BITS, TRI_BITS>(isct, edge, tria);
  return isct;
  }

inline ga::vector<typename fixint::bit_int<161>::type> triangle_edge_intersection_coordinate_fixed(const ga::vector<int32_t>& v0, const ga::vector<int32_t>& v1, const ga::vector<int32_t>& v2, const ga::vector<int32_t>& e0, const ga::vector<int32_t>& e1)
  {
  ga::vector<fixint::limb_int<1>> V0, V1, V2, E0, E1;
  V0.e0 = fixint::limb_int<1>(v0.e0);
  V0.e1 = fixint::limb_int<1>(v0.e1);
  V0.e2 = fixint::limb_int<1>(v0.e2);
  V0.e3 = fixint::limb_int<1>(v0.e3);
  V1.e0 = fixint::limb_int<1>(v1.e0);
  V1.e1 = fixint::limb_int<1>(v1.e1);
  V1.e2 = fixint::limb_int<1>(v1.e2);
  V1.e3 = fixint::limb_int<1>(v1.e3);
  V2.e0 = fixint::limb_int<1>(v2.e0);
  V2.e1 = fixint::limb_int<1>(v2.e1);
  V2.e2 = fixint::limb_int<1>(v2.e2);
  V2.e3 = fixint::limb_int<1>(v2.e3);
  E0.e0 = fixint::limb_int<1>(e0.e0);
  E0.e1 = fixint::limb_int<1>(e0.e1);
  E0.e2 = fixint::limb_int<1>(e0.e2);
  E0.e3 = fixint::limb_int<1>(e0.e3);
  E1.e0 = fixint::limb_int<1>(e1.e0);
  E1.e1 = fixint::limb_int<1>(e1.e1);
  E1.e2 = fixint::limb_int<1>(e1.e2);
  E1.e3 = fixint::limb_int<1>(e1.e3);
  return triangle_edge_intersection_coordinate_fixed(V0, V1, V2, E0, E1);
  }


inline ga::vector<typename fixint::bit_int<291>::type> triangle_triangle_triangle_intersection_coordinate_fixed(const ga::vector<fixint::limb_int<1>>& u0, const ga::vector<fixint::limb_int<1>>& u1, const ga::vector<fixint::limb_int<1>>& u2, const ga::vector<fixint::limb_int<1>>& v0, const ga::vector<fixint::limb_int<1>>& v1, const ga::vector<fixint::limb_int<1>>& v2, const ga::vector<fixint::limb_int<1>>& w0, const ga::vector<fixint::limb_int<1>>& w1, const ga::vector<fixint::limb_int<1>>& w2)
  {
  const static int IN_BITS = 31;
  const static int EXT2_UP_BITS = 2 * IN_BITS + 1; //63
  const static int EXT3_UP_BITS = EXT2_UP_BITS + IN_BITS + 2; //96
  const static int EXT2_DN_BITS = 2 * EXT3_UP_BITS + 1; //193
  const static int ISCT_BITS = EXT2_DN_BITS + EXT3_UP_BITS + 2; //291

  ga::bivector<typename fixint::bit_int<EXT2_UP_BITS>::type> temp;
  ga::trivector<typename fixint::bit_int<EXT3_UP_BITS>::type> tria[3];
  join<IN_BITS, IN_BITS>(temp, u0, u1);
  join<EXT2_UP_BITS, IN_BITS>(tria[0], temp, u2);
  join<IN_BITS, IN_BITS>(temp, v0, v1);
  join<EXT2_UP_BITS, IN_BITS>(tria[1], temp, v2);
  join<IN_BITS, IN_BITS>(temp, w0, w1);
  join<EXT2_UP_BITS, IN_BITS>(tria[2], temp, w2);

  ga::bivector<typename fixint::bit_int<EXT2_DN_BITS>::type> isct_edge;
  ga::vector<typename fixint::bit_int<ISCT_BITS>::type> isct;
  meet<EXT3_UP_BITS, EXT3_UP_BITS>(isct_edge, tria[0], tria[1]);
  meet<EXT2_DN_BITS, EXT3_UP_BITS>(isct, isct_edge, tria[2]);
  return isct;
  }


inline ga::vector<typename fixint::bit_int<579>::type> triangle_triangle_triangle_intersection_coordinate_fixed(const ga::vector<fixint::limb_int<2>>& u0, const ga::vector<fixint::limb_int<2>>& u1, const ga::vector<fixint::limb_int<2>>& u2, const ga::vector<fixint::limb_int<2>>& v0, const ga::vector<fixint::limb_int<2>>& v1, const ga::vector<fixint::limb_int<2>>& v2, const ga::vector<fixint::limb_int<2>>& w0, const ga::vector<fixint::limb_int<2>>& w1, const ga::vector<fixint::limb_int<2>>& w2)
  {
  const static int IN_BITS = 63;
  const static int EXT2_UP_BITS = 2 * IN_BITS + 1; //127
  const static int EXT3_UP_BITS = EXT2_UP_BITS + IN_BITS + 2; //192
  const static int EXT2_DN_BITS = 2 * EXT3_UP_BITS + 1; //385
  const static int ISCT_BITS = EXT2_DN_BITS + EXT3_UP_BITS + 2; //579

  ga::bivector<typename fixint::bit_int<EXT2_UP_BITS>::type> temp;
  ga::trivector<typename fixint::bit_int<EXT3_UP_BITS>::type> tria[3];
  join<IN_BITS, IN_BITS>(temp, u0, u1);
  join<EXT2_UP_BITS, IN_BITS>(tria[0], temp, u2);
  join<IN_BITS, IN_BITS>(temp, v0, v1);
  join<EXT2_UP_BITS, IN_BITS>(tria[1], temp, v2);
  join<IN_BITS, IN_BITS>(temp, w0, w1);
  join<EXT2_UP_BITS, IN_BITS>(tria[2], temp, w2);

  ga::bivector<typename fixint::bit_int<EXT2_DN_BITS>::type> isct_edge;
  ga::vector<typename fixint::bit_int<ISCT_BITS>::type> isct;
  meet<EXT3_UP_BITS, EXT3_UP_BITS>(isct_edge, tria[0], tria[1]);
  meet<EXT2_DN_BITS, EXT3_UP_BITS>(isct, isct_edge, tria[2]);
  return isct;
  }


inline ga::vector<typename fixint::bit_int<291>::type> triangle_triangle_triangle_intersection_coordinate_fixed(const ga::vector<int32_t>& u0, const ga::vector<int32_t>& u1, const ga::vector<int32_t>& u2, const ga::vector<int32_t>& v0, const ga::vector<int32_t>& v1, const ga::vector<int32_t>& v2, const ga::vector<int32_t>& w0, const ga::vector<int32_t>& w1, const ga::vector<int32_t>& w2)
  {
  ga::vector<fixint::limb_int<1>> U0, U1, U2, V0, V1, V2, W0, W1, W2;
  U0.e0 = fixint::limb_int<1>(u0.e0);
  U0.e1 = fixint::limb_int<1>(u0.e1);
  U0.e2 = fixint::limb_int<1>(u0.e2);
  U0.e3 = fixint::limb_int<1>(u0.e3);
  U1.e0 = fixint::limb_int<1>(u1.e0);
  U1.e1 = fixint::limb_int<1>(u1.e1);
  U1.e2 = fixint::limb_int<1>(u1.e2);
  U1.e3 = fixint::limb_int<1>(u1.e3);
  U2.e0 = fixint::limb_int<1>(u2.e0);
  U2.e1 = fixint::limb_int<1>(u2.e1);
  U2.e2 = fixint::limb_int<1>(u2.e2);
  U2.e3 = fixint::limb_int<1>(u2.e3);
  V0.e0 = fixint::limb_int<1>(v0.e0);
  V0.e1 = fixint::limb_int<1>(v0.e1);
  V0.e2 = fixint::limb_int<1>(v0.e2);
  V0.e3 = fixint::limb_int<1>(v0.e3);
  V1.e0 = fixint::limb_int<1>(v1.e0);
  V1.e1 = fixint::limb_int<1>(v1.e1);
  V1.e2 = fixint::limb_int<1>(v1.e2);
  V1.e3 = fixint::limb_int<1>(v1.e3);
  V2.e0 = fixint::limb_int<1>(v2.e0);
  V2.e1 = fixint::limb_int<1>(v2.e1);
  V2.e2 = fixint::limb_int<1>(v2.e2);
  V2.e3 = fixint::limb_int<1>(v2.e3);
  W0.e0 = fixint::limb_int<1>(w0.e0);
  W0.e1 = fixint::limb_int<1>(w0.e1);
  W0.e2 = fixint::limb_int<1>(w0.e2);
  W0.e3 = fixint::limb_int<1>(w0.e3);
  W1.e0 = fixint::limb_int<1>(w1.e0);
  W1.e1 = fixint::limb_int<1>(w1.e1);
  W1.e2 = fixint::limb_int<1>(w1.e2);
  W1.e3 = fixint::limb_int<1>(w1.e3);
  W2.e0 = fixint::limb_int<1>(w2.e0);
  W2.e1 = fixint::limb_int<1>(w2.e1);
  W2.e2 = fixint::limb_int<1>(w2.e2);
  W2.e3 = fixint::limb_int<1>(w2.e3);
  return triangle_triangle_triangle_intersection_coordinate_fixed(U0, U1, U2, V0, V1, V2, W0, W1, W2);
  }
