#pragma once

#include <algorithm>
#include <stdint.h>
#include <vector>
#include <string>
#include <cassert>
#include <iostream>

namespace fixint
  {

  template <int N>
  class limb_int
    {
    public:
      uint32_t limbs[N];

      limb_int(int32_t i);
      limb_int(int64_t i);
      limb_int();
    };

  template <int N>
  inline bool is_negative(const limb_int<N>& a);

  template<int Nout, int Nin>
  inline void promote(limb_int<Nout>& out, const limb_int<Nin>& in);

  template<int Nout, int Nlhs, int Nrhs>
  inline void add(limb_int<Nout>& out, limb_int<Nlhs> lhs, limb_int<Nrhs> rhs);

  template <int N>
  inline limb_int<N> neg(const limb_int<N>& a);

  template<int Nout, int Nlhs, int Nrhs>
  inline void sub(limb_int<Nout>& out, limb_int<Nlhs> lhs, limb_int<Nrhs> rhs);

  template <int N>
  inline void div(limb_int<N>& result, int32_t& remainder, limb_int<N> a, int32_t c);

  template <int N>
  inline bool is_zero(const limb_int<N>& a);

  template <int N>
  inline std::string to_string(limb_int<N> a);

  template <int Nout, int Nlhs, int Nrhs>
  inline void mul(limb_int<Nout>& out, limb_int<Nlhs> lhs, limb_int<Nrhs> rhs);

  template <int Nout, int Nin>
  inline void mul(limb_int<Nout>& out, limb_int<Nin> in, int32_t s);

  template <int N>
  inline void from_string(limb_int<N>& out, const char* s);

  template <int N>
  double to_double(limb_int<N> a);

  template <typename T, int N>
  T to_(limb_int<N> a);

  template <int N>
  inline std::ostream& operator << (std::ostream& out, const limb_int<N>& a);

  template <int Nbits>
  struct bit_int
    {
    typedef limb_int<((Nbits - 1) / 32) + 1> type;
    };

  namespace details
    {
    const static uint32_t limb_bit_size = 32;
    const static uint32_t sign_bit_offset = limb_bit_size - 1;
    const static uint32_t mask_sign_bit = uint32_t(1) << sign_bit_offset;
    const static uint32_t zero_mask_pattern = 0;
    const static uint32_t one_mask_pattern = uint32_t(-1);

    template <int N>
    inline uint32_t _sign_limb(const limb_int<N>& a, int n)
      {
      return (a.limbs[n] & mask_sign_bit) ? one_mask_pattern : zero_mask_pattern;
      }

    inline uint32_t _add_two_limbs(uint32_t& result, uint32_t a, uint32_t b, uint32_t carry)
      {
      uint32_t temp;
      if (carry == 0)
        {
        temp = a + b;
        if (temp < a)
          carry = 1;
        }
      else
        {
        carry = 1;
        temp = a + b + carry;
        if (temp > a)
          carry = 0;
        }
      result = temp;
      return carry;
      }


    inline void _div_two_limbs(uint32_t& result, uint32_t& rest, uint32_t hi, uint32_t lo, uint32_t c)
      {
      union
        {
        struct
          {
          uint32_t low;  // 32 bits
          uint32_t high; // 32 bits
          } u_;

        uint64_t u;       // 64 bits
        } hilo;

      hilo.u_.high = hi;
      hilo.u_.low = lo;

      result = uint32_t(hilo.u / c);
      rest = uint32_t(hilo.u % c);
      }


    inline int _find_leading_bit_in_limb_unsigned(uint32_t x)
      {
      if (x == 0)
        return -1;
      uint32_t bit = limb_bit_size - 1;
      while ((x & mask_sign_bit) == 0)
        {
        x = x << 1;
        --bit;
        }
      return bit;
      }

    /*
    this method looks for the highest set bit

    result:
    if 'a' is not zero:
    return value - true
    'limb_id'    - the index of a limb <0..N-1>
    'index'      - the index of this set bit in the word <0..limb_bit_size)

    if 'a' is zero:
    return value - false
    both 'limb_id' and 'index' are zero
    */
    template <int N>
    inline bool _find_leading_bit_unsigned(int& limb_id, int& index, const limb_int<N>& a)
      {
      for (limb_id = N - 1; limb_id != 0 && a.limbs[limb_id] == 0; --limb_id);
      if (limb_id == 0 && a.limbs[limb_id] == 0)
        {
        index = 0;
        return false;
        }
      index = _find_leading_bit_in_limb_unsigned(a.limbs[limb_id]);
      return true;
      }


    /*
    multiplication: result_high:result_low = a*b
    */
    inline void _mul_two_limbs(uint32_t& result_high, uint32_t& result_low, uint32_t a, uint32_t b)
      {
      union uint_
        {
        struct
          {
          uint32_t low;
          uint32_t high;
          } u_;
        uint64_t u;
        } res;
      res.u = uint64_t(a)*uint64_t(b);
      result_high = res.u_.high;
      result_low = res.u_.low;
      }

    /*!
    this method adds only two unsigned 32-bit limbs to the existing value
    and these limbs begin on the 'index' position
    (it's used in the multiplication algorithm)

    index should be equal or smaller than N-2 (index <= N-2)
    lo - lower limb, hi - higher limb

    for example if we've got N equal 4 and:
    limbs[0] = 3
    limbs[1] = 4
    limbs[2] = 5
    limbs[3] = 6
    then let
    lo = 10
    hi = 20
    and
    index = 1

    the result of this method will be:
    limbs[0] = 3
    limbs[1] = 4 + lo = 14
    limbs[2] = 5 + hi = 25
    limbs[3] = 6

    and no carry at the end of limbs[3]

    (of course if there was a carry in limbs[2] (5+20) then
    this carry would be passed to the limbs[3] etc.)
    */
    inline uint32_t _add_two_ints(uint32_t* limbs, uint32_t N, uint32_t hi, uint32_t lo, int index)
      {
      assert(((uint32_t)index + 1) < N);
      uint32_t carry;
      carry = _add_two_limbs(limbs[index], lo, limbs[index], 0);
      carry = _add_two_limbs(limbs[index + 1], hi, limbs[index + 1], carry);

      for (int i = index + 2; i < (int)N && carry; ++i)
        carry = _add_two_limbs(limbs[i], 0, limbs[i], carry);

      return carry;
      }

    /*!
    this method adds one limb (at a specific position)
    and returns a carry (if it was there)

    if we've got (N=3):
    limbs[0] = 10;
    limbs[1] = 30;
    limbs[2] = 5;
    and we call:
    _add_int(2,1)
    then it'll be:
    limbs[0] = 10;
    limbs[1] = 30 + 2;
    limbs[2] = 5;

    if there was a carry from limbs[2] it would be returned
    */
    inline uint32_t _add_int(uint32_t* limbs, uint32_t N, uint32_t value, int index)
      {
      assert(((uint32_t)index) < N);
      uint32_t carry = _add_two_limbs(limbs[index], value, limbs[index], 0);

      for (int i = index + 1; i < (int)N && carry; ++i)
        carry = _add_two_limbs(limbs[i], 0, limbs[i], carry);

      return carry;
      }

    inline void _mul3(uint32_t* result, uint32_t sz, const uint32_t* ss1, const uint32_t* ss2, uint32_t x1start, uint32_t x1size, uint32_t x2start, uint32_t x2size)
      {
      uint32_t hi, lo;
      for (uint32_t i = 0; i < sz; ++i)
        result[i] = 0;
      if (x1size == 0 || x2size == 0)
        return;
      for (uint32_t x1 = x1start; x1 < x1size; ++x1)
        {
        for (uint32_t x2 = x2start; x2 < x2size; ++x2)
          {
          _mul_two_limbs(hi, lo, ss1[x1], ss2[x2]);
          if (hi != 0)
            {
            auto ind = x2 + x1;
            if (sz > ind + 1)
              _add_two_ints(result, sz, hi, lo, ind);
            else
              {
              _add_int(result, sz, lo, ind);
              }
            }
          else
            {
            if (lo != 0)
              {
              auto ind = x2 + x1;
              auto carry = _add_int(result, sz, lo, ind);
              if (carry)
                {
                _add_int(result, sz, carry, ind + 1);
                }
              }
            }
          }
        }
      }

    inline void _mul2(uint32_t* result, uint32_t sz, const uint32_t* ss1, const uint32_t* ss2, uint32_t sz1, uint32_t sz2)
      {
      uint32_t x1start = 0;
      uint32_t x1size = uint32_t(sz1);
      uint32_t x2start = 0;
      uint32_t x2size = uint32_t(sz2);

      for (x1size = uint32_t(sz1); x1size > 0 && ss1[x1size - 1] == 0; --x1size);
      for (x2size = uint32_t(sz2); x2size > 0 && ss2[x2size - 1] == 0; --x2size);
      for (x1start = 0; x1start < x1size && ss1[x1start] == 0; ++x1start);
      for (x2start = 0; x2start < x2size && ss2[x2start] == 0; ++x2start);

      _mul3(result, sz, ss1, ss2, x1start, x1size, x2start, x2size);
      }
    }


  template <int N>
  inline limb_int<N>::limb_int(int32_t i)
    {
    limbs[0] = uint32_t(i);
    if constexpr (N > 1)
      {
      auto fill = (limbs[0] & details::mask_sign_bit) ? details::one_mask_pattern : details::zero_mask_pattern;
      for (int j = 1; j < N; ++j)
        limbs[j] = fill;
      }
    }

  template <int N>
  inline limb_int<N>::limb_int(int64_t i)
    {
    union uint_
      {
      struct
        {
        uint32_t low;
        uint32_t high;
        } u_;
      uint64_t u;
      } res;
    res.u = uint64_t(i);
    limbs[0] = res.u_.low;
    if constexpr (N > 1)
      {
      limbs[1] = res.u_.high;
      auto fill = (limbs[1] & details::mask_sign_bit) ? details::one_mask_pattern : details::zero_mask_pattern;
      for (int j = 2; j < N; ++j)
        limbs[j] = fill;
      }
    }

  template <int N>
  inline limb_int<N>::limb_int()
    {
    }

  template <int N>
  inline bool is_negative(const limb_int<N>& a)
    {
    return details::_sign_limb(a, N - 1) != details::zero_mask_pattern;
    }

  template<int Nout, int Nin>
  inline void promote(limb_int<Nout>& out, const limb_int<Nin>& in)
    {
    assert(Nout >= Nin);
    for (int i = 0; i < Nin; ++i)
      out.limbs[i] = in.limbs[i];

    if constexpr (Nout > Nin)
      {
      auto fill = details::_sign_limb(out, Nin - 1);
      for (int i = Nin; i < Nout; ++i)
        out.limbs[i] = fill;
      }
    }


  template<int Nout, int Nlhs, int Nrhs>
  inline void add(limb_int<Nout>& out, limb_int<Nlhs> lhs, limb_int<Nrhs> rhs)
    {
    assert(Nout >= Nlhs && Nout >= Nrhs);
    auto pattern_lhs = details::_sign_limb(lhs, Nlhs - 1);
    auto pattern_rhs = details::_sign_limb(rhs, Nrhs - 1);
    uint32_t carry = 0;
    int Nmax;
    if constexpr (Nlhs == Nrhs)
      {
      Nmax = Nlhs;
      for (int k = 0; k < Nrhs; ++k)
        carry = details::_add_two_limbs(out.limbs[k], lhs.limbs[k], rhs.limbs[k], carry);
      }
    else if constexpr (Nlhs > Nrhs)
      {
      Nmax = Nlhs;
      for (int k = 0; k < Nrhs; ++k)
        carry = details::_add_two_limbs(out.limbs[k], lhs.limbs[k], rhs.limbs[k], carry);
      for (int k = Nrhs; k < Nlhs; ++k)
        carry = details::_add_two_limbs(out.limbs[k], lhs.limbs[k], pattern_rhs, carry);
      }
    else
      {
      Nmax = Nrhs;
      for (int k = 0; k < Nlhs; ++k)
        carry = details::_add_two_limbs(out.limbs[k], lhs.limbs[k], rhs.limbs[k], carry);
      for (int k = Nlhs; k < Nrhs; ++k)
        carry = details::_add_two_limbs(out.limbs[k], pattern_lhs, rhs.limbs[k], carry);
      }
    if (Nout > Nmax)
      {
      if (!pattern_rhs && !pattern_lhs)
        {
        for (int i = Nmax; i < Nout; ++i)
          out.limbs[i] = details::zero_mask_pattern;
        }
      else if (pattern_rhs && pattern_lhs)
        {
        for (int i = Nmax; i < Nout; ++i)
          out.limbs[i] = details::one_mask_pattern;
        }
      else
        {
        auto sign = details::_sign_limb(out, Nmax - 1);
        auto fill = (sign) ? details::one_mask_pattern : details::zero_mask_pattern;
        for (int i = Nmax; i < Nout; ++i)
          out.limbs[i] = fill;
        }
      }
    }

  template <int N>
  inline limb_int<N> neg(const limb_int<N>& a)
    {
    limb_int<N> n(a);
    for (int i = 0; i < N; ++i)
      n.limbs[i] ^= details::one_mask_pattern;
    add(n, n, limb_int<N>(1));
    return n;
    }

  template<int Nout, int Nlhs, int Nrhs>
  inline void sub(limb_int<Nout>& out, limb_int<Nlhs> lhs, limb_int<Nrhs> rhs)
    {
    rhs = neg(rhs);
    add(out, lhs, rhs);
    }


  template <int N>
  inline void div(limb_int<N>& result, int32_t& remainder, limb_int<N> a, int32_t c)
    {
    if (c == 1)
      {
      remainder = 0;
      result = a;
      return;
      }
    bool negative = is_negative(a);
    bool neg_div = false;
    if (c < 0)
      {
      neg_div = true;
      c = -c;
      }
    if (negative)
      a = neg(a);
    uint32_t rest = 0;
    int i;
    for (i = N - 1; i > 0 && a.limbs[i] == 0; --i)
      result.limbs[i] = 0;
    for (; i >= 0; --i)
      details::_div_two_limbs(result.limbs[i], rest, rest, a.limbs[i], c);
    if (negative != neg_div)
      {
      if (!is_negative(result))
        result = neg(result);
      }
    if (negative)
      remainder = -int32_t(rest);
    else
      remainder = int32_t(rest);
    }


  template <int N>
  inline bool is_zero(const limb_int<N>& a)
    {
    for (int i = 0; i < N; ++i)
      {
      if (a.limbs[i])
        return false;
      }
    return true;
    }

  template <int N>
  inline std::string to_string(limb_int<N> a)
    {
    bool negative = is_negative(a);
    std::string result;
    int limb_id, index;
    if (!details::_find_leading_bit_unsigned(limb_id, index, a))
      {
      result = '0';
      return result;
      }
    double digits_d = limb_id;
    digits_d *= details::limb_bit_size;
    digits_d += index + 1;
    digits_d *= 0.301029995663981195;
    uint32_t digits = static_cast<uint32_t>(digits_d) + 3;
    if (result.capacity() < digits)
      result.reserve(digits);
    int32_t rest;
    char ch;
    do
      {
      div(a, rest, a, 10);
      ch = static_cast<char>(std::abs(rest) + '0');
      result.insert(result.end(), ch);
      } while (!is_zero(a));
      if (negative)
        result.insert(result.end(), '-');
      std::reverse(result.begin(), result.end());
      return result;
    }


  template <int Nout, int Nlhs, int Nrhs>
  inline void mul(limb_int<Nout>& out, limb_int<Nlhs> lhs, limb_int<Nrhs> rhs)
    {
    assert(Nout >= Nlhs);
    assert(Nout >= Nrhs);
    bool neg_lhs = is_negative(lhs);
    bool neg_rhs = is_negative(rhs);
    if (neg_lhs)
      lhs = neg(lhs);
    if (neg_rhs)
      rhs = neg(rhs);
    details::_mul2(out.limbs, Nout, lhs.limbs, rhs.limbs, Nlhs, Nrhs);
    if (neg_lhs != neg_rhs)
      {
      if (!is_negative(out))
        out = neg(out);
      }
    }

  template <int Nout, int Nin>
  inline void mul(limb_int<Nout>& out, limb_int<Nin> in, int32_t s)
    {    
    mul(out, in, limb_int<1>(s));
    }

  template <int N>
  inline void from_string(limb_int<N>& out, const char* s)
    {
    bool negative = false;
    out = limb_int<N>(0);
    while ((*s == ' ') || (*s == '\t') || (*s == 13) || (*s == '\n'))
      ++s;
    if (*s == '-')
      {
      negative = true;
      ++s;
      while ((*s == ' ') || (*s == '\t') || (*s == 13) || (*s == '\n'))
        ++s;
      }    
    for (; (*s >= '0' && *s <= '9'); ++s)
      {
      int32_t d = *s - '0';
      mul(out, out, 10);
      if (negative)
        sub(out, out, limb_int<1>(d));
      else
        add(out, out, limb_int<1>(d));
      }
    }

  template <int N>
  double to_double(limb_int<N> a)
    {
    const static uint64_t factor = uint64_t(1) << uint64_t(details::limb_bit_size);
    bool negative = is_negative(a);
    if (negative)
      a = neg(a);
    double d = double(a.limbs[N - 1]);
    for (int k = N - 2; k >= 0; --k)
      {
      d *= double(factor);
      d += double(a.limbs[k]);
      }
    return negative ? -d : d;
    }

  template <typename T, int N>
  T to_(limb_int<N> a)
    {
    const static uint64_t factor = uint64_t(1) << uint64_t(details::limb_bit_size);
    bool negative = is_negative(a);
    if (negative)
      a = neg(a);
    T d = T(a.limbs[N - 1]);
    for (int k = N - 2; k >= 0; --k)
      {
      d *= T(factor);
      d += T(a.limbs[k]);
      }
    return negative ? -d : d;
    }

  template <int N>
  inline std::ostream& operator << (std::ostream& out, const limb_int<N>& a)
    {
    if ((out.flags() & out.basefield) == out.hex)
      {
      out << "[" << std::hex << a.limbs[0];
      for (size_t k = 1; k < N; ++k)
        out << ";" << std::hex << a.limbs[k];
      out << "]";
      }
    else
      out << to_string(a);
    return out;
    }

  }