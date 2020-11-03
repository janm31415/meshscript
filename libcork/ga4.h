#pragma once

/*
A Survey of Geometric Algebra and Geometric Calculus. Alan Macdonald.
http://faculty.luther.edu/~macdonal/GA&GC.pdf

This code implements G4, a geometric algebra for the homogeneous model of 3D space (so it is a 4D algebra).

*/

#include <cmath>

namespace ga
  {

  template <class TType>
  struct vector
    {
    TType e0;
    TType e1;
    TType e2;
    TType e3;
    };

  template <class TType>
  struct bivector
    {
    TType e0e1;
    TType e0e2;
    TType e0e3;
    TType e1e2;
    TType e1e3;
    TType e2e3;
    };

  template <class TType>
  struct trivector
    {
    TType e0e1e2;
    TType e0e1e3;
    TType e0e2e3;
    TType e1e2e3;
    };

  template <class TType>
  struct vector4
    {
    TType e0e1e2e3;
    };

  template <class TType>
  vector<TType> operator - (const vector<TType>& v)
    {
    vector<TType> res;
    res.e0 = -v.e0;
    res.e1 = -v.e1;
    res.e2 = -v.e2;
    res.e3 = -v.e3;
    return res;
    }

  template <class TType>
  bivector<TType> operator - (const bivector<TType>& v)
    {
    bivector<TType> res;
    res.e0e1 = -v.e0e1;
    res.e0e2 = -v.e0e2;
    res.e0e3 = -v.e0e3;
    res.e1e2 = -v.e1e2;
    res.e1e3 = -v.e1e3;
    res.e2e3 = -v.e2e3;
    return res;
    }

  template <class TType>
  trivector<TType> operator - (const trivector<TType>& v)
    {
    trivector<TType> res;
    res.e0e1e2 = -v.e0e1e2;
    res.e0e1e3 = -v.e0e1e3;
    res.e0e2e3 = -v.e0e2e3;
    res.e1e2e3 = -v.e1e2e3;
    return res;
    }

  template <class TType>
  vector4<TType> operator - (const vector4<TType>& v)
    {
    vector4<TType> res;
    res.e0e1e2e3 = -v.e0e1e2e3;    
    return res;
    }

  template <class TType>
  vector<TType> abs(const vector<TType>& v)
    {
    vector<TType> res;
    res.e0 = fabs(v.e0);
    res.e1 = fabs(v.e1);
    res.e2 = fabs(v.e2);
    res.e3 = fabs(v.e3);
    return res;
    }

  template <class TType>
  bivector<TType> abs(const bivector<TType>& v)
    {
    bivector<TType> res;
    res.e0e1 = fabs(v.e0e1);
    res.e0e2 = fabs(v.e0e2);
    res.e0e3 = fabs(v.e0e3);
    res.e1e2 = fabs(v.e1e2);
    res.e1e3 = fabs(v.e1e3);
    res.e2e3 = fabs(v.e2e3);
    return res;
    }

  template <class TType>
  trivector<TType> abs(const trivector<TType>& v)
    {
    trivector<TType> res;
    res.e0e1e2 = fabs(v.e0e1e2);
    res.e0e1e3 = fabs(v.e0e1e3);
    res.e0e2e3 = fabs(v.e0e2e3);
    res.e1e2e3 = fabs(v.e1e2e3);
    return res;
    }

  template <class TType>
  vector4<TType> abs(const vector4<TType>& v)
    {
    vector4<TType> res;
    res.e0e1e2e3 = fabs(v.e0e1e2e3);
    return res;
    }

  template <class TType>
  vector<TType> operator * (const vector<TType>& v, const TType& s)
    {
    vector<TType> res;
    res.e0 = s*v.e0;
    res.e1 = s*v.e1;
    res.e2 = s*v.e2;
    res.e3 = s*v.e3;
    return res;
    }

  template <class TType>
  bivector<TType> operator * (const bivector<TType>& v, const TType& s)
    {
    bivector<TType> res;
    res.e0e1 = s*v.e0e1;
    res.e0e2 = s*v.e0e2;
    res.e0e3 = s*v.e0e3;
    res.e1e2 = s*v.e1e2;
    res.e1e3 = s*v.e1e3;
    res.e2e3 = s*v.e2e3;
    return res;
    }

  template <class TType>
  trivector<TType> operator * (const trivector<TType>& v, const TType& s)
    {
    trivector<TType> res;
    res.e0e1e2 = s*v.e0e1e2;
    res.e0e1e3 = s*v.e0e1e3;
    res.e0e2e3 = s*v.e0e2e3;
    res.e1e2e3 = s*v.e1e2e3;
    return res;
    }

  template <class TType>
  vector4<TType> operator * (const vector4<TType>& v, const TType& s)
    {
    vector4<TType> res;
    res.e0e1e2e3 = s*v.e0e1e2e3;
    return res;
    }

  template <class TType>
  vector<TType> operator * (const TType& s, const vector<TType>& v)
    {
    vector<TType> res;
    res.e0 = s*v.e0;
    res.e1 = s*v.e1;
    res.e2 = s*v.e2;
    res.e3 = s*v.e3;
    return res;
    }

  template <class TType>
  bivector<TType> operator * (const TType& s, const bivector<TType>& v)
    {
    bivector<TType> res;
    res.e0e1 = s*v.e0e1;
    res.e0e2 = s*v.e0e2;
    res.e0e3 = s*v.e0e3;
    res.e1e2 = s*v.e1e2;
    res.e1e3 = s*v.e1e3;
    res.e2e3 = s*v.e2e3;
    return res;
    }

  template <class TType>
  trivector<TType> operator * (const TType& s, const trivector<TType>& v)
    {
    trivector<TType> res;
    res.e0e1e2 = s*v.e0e1e2;
    res.e0e1e3 = s*v.e0e1e3;
    res.e0e2e3 = s*v.e0e2e3;
    res.e1e2e3 = s*v.e1e2e3;
    return res;
    }

  template <class TType>
  vector4<TType> operator * (const TType& s, const vector4<TType>& v)
    {
    vector4<TType> res;
    res.e0e1e2e3 = s*v.e0e1e2e3;
    return res;
    }


  template <class TType>
  vector<TType> operator / (const vector<TType>& v, const TType& s)
    {
    vector<TType> res;
    res.e0 = v.e0/s;
    res.e1 = v.e1/s;
    res.e2 = v.e2/s;
    res.e3 = v.e3/s;
    return res;
    }

  template <class TType>
  bivector<TType> operator / (const bivector<TType>& v, const TType& s)
    {
    bivector<TType> res;
    res.e0e1 = v.e0e1/s;
    res.e0e2 = v.e0e2/s;
    res.e0e3 = v.e0e3/s;
    res.e1e2 = v.e1e2/s;
    res.e1e3 = v.e1e3/s;
    res.e2e3 = v.e2e3/s;
    return res;
    }

  template <class TType>
  trivector<TType> operator / (const trivector<TType>& v, const TType& s)
    {
    trivector<TType> res;
    res.e0e1e2 = v.e0e1e2/s;
    res.e0e1e3 = v.e0e1e3/s;
    res.e0e2e3 = v.e0e2e3/s;
    res.e1e2e3 = v.e1e2e3/s;
    return res;
    }

  template <class TType>
  vector4<TType> operator / (const vector4<TType>& v, const TType& s)
    {
    vector4<TType> res;
    res.e0e1e2e3 = v.e0e1e2e3/s;
    return res;
    }

  template <class TType>
  vector<TType> operator + (const vector<TType>& v1, const vector<TType>& v2)
    {
    vector<TType> res(v1);
    res.e0 += v2.e0;
    res.e1 += v2.e1;
    res.e2 += v2.e2;
    res.e3 += v2.e3;
    return res;
    }

  template <class TType>
  bivector<TType> operator + (const bivector<TType>& v1, const bivector<TType>& v2)
    {
    bivector<TType> res(v1);
    res.e0e1 += v2.e0e1;
    res.e0e2 += v2.e0e2;
    res.e0e3 += v2.e0e3;
    res.e1e2 += v2.e1e2;
    res.e1e3 += v2.e1e3;
    res.e2e3 += v2.e2e3;
    return res;
    }

  template <class TType>
  trivector<TType> operator + (const trivector<TType>& v1, const trivector<TType>& v2)
    {
    trivector<TType> res(v1);
    res.e0e1e2 += v2.e0e1e2;
    res.e0e1e3 += v2.e0e1e3;
    res.e0e2e3 += v2.e0e2e3;
    res.e1e2e3 += v2.e1e2e3;
    return res;
    }

  template <class TType>
  vector4<TType> operator + (const vector4<TType>& v1, const vector4<TType>& v2)
    {
    vector4<TType> res(v1);
    res.e0e1e2e3 += v2.e0e1e2e3;
    return res;
    }


  template <class TType>
  vector<TType> operator - (const vector<TType>& v1, const vector<TType>& v2)
    {
    vector<TType> res(v1);
    res.e0 -= v2.e0;
    res.e1 -= v2.e1;
    res.e2 -= v2.e2;
    res.e3 -= v2.e3;
    return res;
    }

  template <class TType>
  bivector<TType> operator - (const bivector<TType>& v1, const bivector<TType>& v2)
    {
    bivector<TType> res(v1);
    res.e0e1 -= v2.e0e1;
    res.e0e2 -= v2.e0e2;
    res.e0e3 -= v2.e0e3;
    res.e1e2 -= v2.e1e2;
    res.e1e3 -= v2.e1e3;
    res.e2e3 -= v2.e2e3;
    return res;
    }

  template <class TType>
  trivector<TType> operator - (const trivector<TType>& v1, const trivector<TType>& v2)
    {
    trivector<TType> res(v1);
    res.e0e1e2 -= v2.e0e1e2;
    res.e0e1e3 -= v2.e0e1e3;
    res.e0e2e3 -= v2.e0e2e3;
    res.e1e2e3 -= v2.e1e2e3;
    return res;
    }

  template <class TType>
  vector4<TType> operator - (const vector4<TType>& v1, const vector4<TType>& v2)
    {
    vector4<TType> res(v1);
    res.e0e1e2e3 -= v2.e0e1e2e3;
    return res;
    }

  template <class TType>
  bool operator == (const vector<TType>& v1, const vector<TType>& v2)
    {
    return (v1.e0 == v2.e0 && v1.e1 == v2.e1 && v1.e2 == v2.e2 && v1.e3 == v2.e3);
    }

  template <class TType>
  bool operator == (const bivector<TType>& v1, const bivector<TType>& v2)
    {
    return (v1.e0e1 == v2.e0e1 && v1.e0e2 == v2.e0e2 && v1.e0e3 == v2.e0e3 && v1.e1e2 == v2.e1e2 && v1.e1e3 == v2.e1e3 && v1.e2e3 == v2.e2e3);
    }

  template <class TType>
  bool operator == (const trivector<TType>& v1, const trivector<TType>& v2)
    {
    return (v1.e0e1e2 == v2.e0e1e2 && v1.e0e1e3 == v2.e0e1e3 && v1.e0e2e3 == v2.e0e2e3 && v1.e1e2e3 == v2.e1e2e3);
    }

  template <class TType>
  bool operator == (const vector4<TType>& v1, const vector4<TType>& v2)
    {
    return (v1.e0e1e2e3 == v2.e0e1e2e3);
    }

  template <class TType>
  bool operator != (const vector<TType>& v1, const vector<TType>& v2)
    {
    return !(v1 == v2);
    }

  template <class TType>
  bool operator != (const bivector<TType>& v1, const bivector<TType>& v2)
    {
    return !(v1 == v2);
    }

  template <class TType>
  bool operator != (const trivector<TType>& v1, const trivector<TType>& v2)
    {
    return !(v1 == v2);
    }

  template <class TType>
  bool operator != (const vector4<TType>& v1, const vector4<TType>& v2)
    {
    return !(v1 == v2);
    }
  /*
  Section 1.5.2/1.5.3 from "A Survey of Geometric Algebra and Geometric Calculus. Alan Macdonald."

  The dual of a multivector A is A* = AI^{-1}. Here I = e0e1e2e3 and I^{-1} = e3e2e1e0
  */

  template <class TType>
  TType dual(const vector4<TType>& v)
    {
    TType res = v.e0e1e2e3;
    return res;
    }

  template <class TType>
  vector<TType> dual(const trivector<TType>& v)
    {
    vector<TType> res;
    res.e0 = v.e1e2e3;
    res.e1 = -v.e0e2e3;
    res.e2 = v.e0e1e3;
    res.e3 = -v.e0e1e2;
    return res;
    }

  template <class TType>
  bivector<TType> dual(const bivector<TType>& v)
    {
    bivector<TType> res;
    res.e0e1 = v.e2e3;
    res.e0e2 = -v.e1e3;
    res.e0e3 = v.e1e2;
    res.e1e2 = v.e0e3;
    res.e1e3 = -v.e0e2;
    res.e2e3 = v.e0e1;
    return res;
    }

  template <class TType>
  trivector<TType> dual(const vector<TType>& v)
    {
    trivector<TType> res;
    res.e0e1e2 = v.e3;
    res.e0e1e3 = -v.e2;
    res.e0e2e3 = v.e1;
    res.e1e2e3 = -v.e0;
    return res;
    }

  template <class TType>
  vector4<TType> dual(const TType& s)
    {
    vector4<TType> res;
    res.e0e1e2e3 = s;
    return res;
    }

  template <class TType>
  TType undo_dual(const vector4<TType>& v)
    {
    TType res = v.e0e1e2e3;
    return res;
    }

  template <class TType>
  vector<TType> undo_dual(const trivector<TType>& v)
    {
    vector<TType> res;
    res.e0 = -v.e1e2e3;
    res.e1 = v.e0e2e3;
    res.e2 = -v.e0e1e3;
    res.e3 = v.e0e1e2;
    return res;
    }

  template <class TType>
  bivector<TType> undo_dual(const bivector<TType>& v)
    {
    bivector<TType> res;
    res.e0e1 = v.e2e3;
    res.e0e2 = -v.e1e3;
    res.e0e3 = v.e1e2;
    res.e1e2 = v.e0e3;
    res.e1e3 = -v.e0e2;
    res.e2e3 = v.e0e1;
    return res;
    }

  template <class TType>
  trivector<TType> undo_dual(const vector<TType>& v)
    {
    trivector<TType> res;
    res.e0e1e2 = -v.e3;
    res.e0e1e3 = v.e2;
    res.e0e2e3 = -v.e1;
    res.e1e2e3 = v.e0;
    return res;
    }

  template <class TType>
  vector4<TType> undo_dual(const TType& s)
    {
    vector4<TType> res;
    res.e0e1e2e3 = s;
    return res;
    }

  /*
  Inner product:
  Let A be a j-vector and B a k-vector. Define A.B = <AB>_{k-j}. Thus if j>k, then A.B = 0
  */

  template <class TType>
  TType inner(const TType& s1, const TType& s2)
    {
    return s1*s2;
    }

  template <class TType>
  TType inner(const vector<TType>& v1, const vector<TType>& v2)
    {
    return v1.e0*v2.e0 + v1.e1*v2.e1 + v1.e2*v2.e2 + v1.e3*v2.e3;
    }

  template <class TType>
  TType inner(const bivector<TType>& v1, const bivector<TType>& v2)
    {
    return v1.e0e1*v2.e0e1 + v1.e0e2*v2.e0e2 + v1.e0e3*v2.e0e3 + v1.e1e2*v2.e1e2 + v1.e1e3*v2.e1e3 + v1.e2e3*v2.e2e3;
    }

  template <class TType>
  TType inner(const trivector<TType>& v1, const trivector<TType>& v2)
    {
    return v1.e0e1e2*v2.e0e1e2 + v1.e0e1e3*v2.e0e1e3 + v1.e0e2e3*v2.e0e2e3 + v1.e1e2e3*v2.e1e2e3;
    }

  template <class TType>
  TType inner(const vector<TType>& v1, const TType& s2)
    {
    return TType(0);
    }

  template <class TType>
  TType inner(const bivector<TType>& v1, const TType& s2)
    {
    return TType(0);
    }

  template <class TType>
  TType inner(const trivector<TType>& v1, const TType& s2)
    {
    return TType(0);
    }

  template <class TType>
  TType inner(const vector4<TType>& v1, const TType& s2)
    {
    return TType(0);
    }

  template <class TType>
  vector<TType> inner(const TType& s1, const vector<TType>& v2)
    {
    return s1*v2;
    }

  template <class TType>
  bivector<TType> inner(const TType& s1, const bivector<TType>& v2)
    {
    return s1*v2;
    }

  template <class TType>
  trivector<TType> inner(const TType& s1, const trivector<TType>& v2)
    {
    return s1*v2;
    }

  template <class TType>
  vector4<TType> inner(const TType& s1, const vector4<TType>& v2)
    {
    return s1*v2;
    }

  template <class TType>
  TType inner(const bivector<TType>& v1, const vector<TType>& v2)
    {
    return TType(0);
    }

  template <class TType>
  TType inner(const trivector<TType>& v1, const vector<TType>& v2)
    {
    return TType(0);
    }

  template <class TType>
  TType inner(const vector4<TType>& v1, const vector<TType>& v2)
    {
    return TType(0);
    }

  template <class TType>
  vector<TType> inner(const vector<TType>& v1, const bivector<TType>& v2)
    {
    vector<TType> res;
    res.e0 = -v1.e1*v2.e0e1 - v1.e2*v2.e0e2 - v1.e3*v2.e0e3;
    res.e1 = v1.e0*v2.e0e1 - v1.e2*v2.e1e2 - v1.e3*v2.e1e3;
    res.e2 = v1.e0*v2.e0e2 + v1.e1*v2.e1e2 - v1.e3*v2.e2e3;
    res.e3 = v1.e0*v2.e0e3 + v1.e1*v2.e1e3 + v1.e2*v2.e2e3;
    return res;
    }

  template <class TType>
  bivector<TType> inner(const vector<TType>& v1, const trivector<TType>& v2)
    {
    bivector<TType> res;
    res.e0e1 = v1.e2*v2.e0e1e2 + v1.e3*v2.e0e1e3;
    res.e0e2 = -v1.e1*v2.e0e1e2 + v1.e3*v2.e0e2e3;
    res.e0e3 = -v1.e1*v2.e0e1e3 - v1.e2*v2.e0e2e3;
    res.e1e2 = v1.e0*v2.e0e1e2 + v1.e3*v2.e1e2e3;
    res.e1e3 = v1.e0*v2.e0e1e3 - v1.e2*v2.e1e2e3;
    res.e2e3 = v1.e0*v2.e0e2e3 + v1.e1*v2.e1e2e3;
    return res;
    }

  template <class TType>
  trivector<TType> inner(const vector<TType>& v1, const vector4<TType>& v2)
    {
    trivector<TType> res;
    res.e0e1e2 = -v1.e3*v2.e0e1e2e3;
    res.e0e1e3 = v1.e2*v2.e0e1e2e3;
    res.e0e2e3 = -v1.e1*v2.e0e1e2e3;
    res.e1e2e3 = v1.e0*v2.e0e1e2e3;
    return res;
    }

  template <class TType>
  TType inner(const trivector<TType>&, const bivector<TType>&)
    {
    return TType(0);
    }

  template <class TType>
  TType inner(const vector4<TType>&, const bivector<TType>&)
    {
    return TType(0);
    }

  template <class TType>
  TType inner(const vector4<TType>&, const trivector<TType>&)
    {
    return TType(0);
    }

  template <class TType>
  vector<TType> inner(const bivector<TType>& v1, const trivector<TType>& v2)
    {
    vector<TType> res;
    res.e0 = -v1.e1e2*v2.e0e1e2 - v1.e1e3*v2.e0e1e3 - v1.e2e3*v2.e0e2e3;
    res.e1 = v1.e0e2*v2.e0e1e2 + v1.e0e3*v2.e0e1e3 - v1.e2e3*v2.e1e2e3;
    res.e2 = -v1.e0e1*v2.e0e1e2 + v1.e0e3*v2.e0e2e3 + v1.e1e3*v2.e1e2e3;
    res.e3 = -v1.e0e1*v2.e0e1e3 - v1.e0e2*v2.e0e2e3 - v1.e1e2*v2.e1e2e3;
    return res;
    }

  template <class TType>
  bivector<TType> inner(const bivector<TType>& v1, const vector4<TType>& v2)
    {
    bivector<TType> res;
    res.e0e1 = -v1.e2e3*v2.e0e1e2e3;
    res.e0e2 = v1.e1e3*v2.e0e1e2e3;
    res.e0e3 = -v1.e1e2*v2.e0e1e2e3;
    res.e1e2 = -v1.e0e3*v2.e0e1e2e3;
    res.e1e3 = v1.e0e2*v2.e0e1e2e3;
    res.e2e3 = -v1.e0e1*v2.e0e1e2e3;
    return res;
    }

  template <class TType>
  vector<TType> inner(const trivector<TType>& v1, const vector4<TType>& v2)
    {
    vector<TType> res;
    return res;
    }

  /*
  Outer product:
  Let A be a j-vector and B a k-vector. Define A^B = <AB>_{j+k}. Thus if j+k>4, then A^B = 0
  */

  template <class TType>
  bivector<TType> outer(const vector<TType>& v1, const vector<TType>& v2)
    {
    bivector<TType> res;
    res.e0e1 = (v1.e0 * v2.e1) - (v2.e0 * v1.e1);
    res.e0e2 = (v1.e0 * v2.e2) - (v2.e0 * v1.e2);
    res.e0e3 = (v1.e0 * v2.e3) - (v2.e0 * v1.e3);
    res.e1e2 = (v1.e1 * v2.e2) - (v2.e1 * v1.e2);
    res.e1e3 = (v1.e1 * v2.e3) - (v2.e1 * v1.e3);
    res.e2e3 = (v1.e2 * v2.e3) - (v2.e2 * v1.e3);
    return res;
    }

  template <class TType>
  trivector<TType> outer(const bivector<TType>& v1, const vector<TType>& v2)
    {
    trivector<TType> res;
    res.e0e1e2 = (v1.e0e1 * v2.e2) - (v1.e0e2 * v2.e1) + (v1.e1e2 * v2.e0);
    res.e0e1e3 = (v1.e0e1 * v2.e3) - (v1.e0e3 * v2.e1) + (v1.e1e3 * v2.e0);
    res.e0e2e3 = (v1.e0e2 * v2.e3) - (v1.e0e3 * v2.e2) + (v1.e2e3 * v2.e0);
    res.e1e2e3 = (v1.e1e2 * v2.e3) - (v1.e1e3 * v2.e2) + (v1.e2e3 * v2.e1);
    return res;
    }

  template <class TType>
  trivector<TType> outer(const vector<TType>& v1, const bivector<TType>& v2)
    {
    trivector<TType> res;
    res.e0e1e2 = (v1.e0 * v2.e1e2) - (v1.e1 * v2.e0e2) + (v1.e2 * v2.e0e1);
    res.e0e1e3 = (v1.e0 * v2.e1e3) - (v1.e1 * v2.e0e3) + (v1.e3 * v2.e0e1);
    res.e0e2e3 = (v1.e0 * v2.e2e3) - (v1.e2 * v2.e0e3) + (v1.e3 * v2.e0e2);
    res.e1e2e3 = (v1.e1 * v2.e2e3) - (v1.e2 * v2.e1e3) + (v1.e3 * v2.e1e2);
    }

  template <class TType>
  vector4<TType> outer(const trivector<TType>& v1, const vector<TType>& v2)
    {
    vector4<TType> res;
    res.e0e1e2e3 = (v1.e0e1e2 * v2.e3) - (v1.e0e1e3 * v2.e2) + (v1.e0e2e3 * v2.e1) - (v1.e1e2e3 * v2.e0);
    return res;
    }

  template <class TType>
  vector4<TType> outer(const bivector<TType>& v1, const bivector<TType>& v2)
    {
    vector4<TType> res;
    res.e0e1e2e3 = (v1.e0e1 * v2.e2e3) - (v1.e0e2 * v2.e1e3) + (v1.e0e3 * v2.e1e2) + (v1.e1e2 * v2.e0e3) - (v1.e1e3 * v2.e0e2) + (v1.e2e3 * v2.e0e1);
    return res;
    }

  template <class TType>
  vector4<TType> outer(const vector<TType>& v1, const trivector<TType>& v2)
    {
    vector4<TType> res;
    res.e0e1e2e3 = (v1.e0 * v2.e1e2e3) - (v1.e1 * v2.e0e2e3) + (v1.e2 * v2.e0e1e3) - (v1.e3 * v2.e0e1e2);
    return res;
    }

  template <class TType>
  vector<TType> outer(const TType& s1, const vector<TType>& v2)
    {
    return v2*s1;
    }

  template <class TType>
  bivector<TType> outer(const TType& s1, const bivector<TType>& v2)
    {
    return v2*s1;
    }

  template <class TType>
  trivector<TType> outer(const TType& s1, const trivector<TType>& v2)
    {
    return v2*s1;
    }

  template <class TType>
  vector4<TType> outer(const TType& s1, const vector4<TType>& v2)
    {
    return v2*s1;
    }

  template <class TType>
  vector<TType> outer(const vector<TType>& v1, const TType& s2)
    {
    return v1*s2;
    }

  template <class TType>
  bivector<TType> outer(const bivector<TType>& v1, const TType& s2)
    {
    return v1*s2;
    }

  template <class TType>
  trivector<TType> outer(const trivector<TType>& v1, const TType& s2)
    {
    return v1*s2;
    }

  template <class TType>
  vector4<TType> outer(const vector4<TType>& v1, const TType& s2)
    {
    return v1*s2;
    }
  }

