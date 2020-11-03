#pragma once

#include <stdint.h>

class quantization
  {
  public:
    quantization(double magnitude, int bits)
      {
      reset(magnitude, bits);
      }

    ~quantization()
      {
      }

    void reset(double magnitude, int bits)
      {
      int max_exponent;
      std::frexp(magnitude, &max_exponent);
      ++max_exponent; // ensure that 2^max_exponent > magnitude
                      // set constants
      _magnify_factor = std::pow((double)2, bits - max_exponent);
      // It is guaranteed that magnitude * _magnify_factor < 2.0^m_bits
      _reshrink_factor = std::pow((double)2, max_exponent - bits);
      }

    double operator()(double number) const
      {
      return _reshrink_factor * (double)((int64_t)(number * _magnify_factor));
      }

    int64_t to_int(double number) const
      {
      return (int64_t)(number * _magnify_factor);
      }

    double to_real(int64_t number) const
      {
      return _reshrink_factor * (double)(number);
      }

    double get_reshrink_factor() const
      {
      return _reshrink_factor;
      }

    double get_magnify_factor() const
      {
      return _magnify_factor;
      }

  private:
    double _magnify_factor;
    double _reshrink_factor;
  };

