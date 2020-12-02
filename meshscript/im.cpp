#include "im.h"

#include <stb_image.h>
#include <stb_image_write.h>

#include <jtk/file_utils.h>

#include <algorithm>
#include <iostream>

bool read_from_file(im& i, const std::string& filename)
  {
  int w, h, nr_of_channels;
  unsigned char* im = stbi_load(filename.c_str(), &w, &h, &nr_of_channels, 4);
  if (im)
    {
    i.texture = jtk::span_to_image(w, h, w, (const uint32_t*)im);
    stbi_image_free(im);
    return true;
    }
  return false;
  }

bool write_to_file(const im& i, const std::string& filename)
  {
  std::string ext = jtk::get_extension(filename);
  if (ext.empty())
    return false;
  std::transform(ext.begin(), ext.end(), ext.begin(), [](char ch) {return (char)::tolower(ch); });
  if (ext == "png")
    {
    if (!stbi_write_png(filename.c_str(), i.texture.width(), i.texture.height(), 4, (void*)i.texture.data(), i.texture.stride() * 4))
      return false;
    }
  else if (ext == "jpg" || ext == "jpeg")
    {
    if (!stbi_write_jpg(filename.c_str(), i.texture.width(), i.texture.height(), 4, (void*)i.texture.data(), 80))
      return false;
    }
  else if (ext == "bmp")
    {
    if (!stbi_write_bmp(filename.c_str(), i.texture.width(), i.texture.height(), 4, (void*)i.texture.data()))
      return false;
    }
  else if (ext == "tga")
    {
    if (!stbi_write_tga(filename.c_str(), i.texture.width(), i.texture.height(), 4, (void*)i.texture.data()))
      return false;
    }
  return true;
  }


void info(const im& i)
  {
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "IMAGE" << std::endl;
  std::cout << "Width: " << i.texture.width() << std::endl;
  std::cout << "Height: " << i.texture.height() << std::endl;
  std::cout << "---------------------------------------" << std::endl;
  }

namespace jtk
  {

  inline void unpack_8bit_to_16bit(const __m128i a, __m128i& b0, __m128i& b1)
    {
    __m128i zero = _mm_setzero_si128();
    b0 = _mm_unpacklo_epi8(a, zero);
    b1 = _mm_unpackhi_epi8(a, zero);
    }

  inline void split_hi_lo(const __m128i a, const __m128i b, __m128i& hi, __m128i& lo)
    {
    hi = _mm_unpackhi_epi8(a, b);
    lo = _mm_unpacklo_epi8(a, b);
    }

  void convolve_col_14641_sse(image<int16_t>& out, const image<uint8_t>& im, bool border = false)
    {
    const int w_chunk = im.stride() / 16;
    out = image<int16_t>(im.width(), im.height());
    __m128i* i0 = (__m128i*)(im.data());
    __m128i* i1 = (__m128i*)(im.data()) + w_chunk * 1;
    __m128i* i2 = (__m128i*)(im.data()) + w_chunk * 2;
    __m128i* i3 = (__m128i*)(im.data()) + w_chunk * 3;
    __m128i* i4 = (__m128i*)(im.data()) + w_chunk * 4;
    __m128i* result_v = (__m128i*)(out.data()) + 4 * w_chunk;
    __m128i* end_input = (__m128i*)(im.data()) + w_chunk * im.height();
    __m128i sixes = _mm_set1_epi16(6);
    __m128i fours = _mm_set1_epi16(4);
    __m128i sixteen = _mm_set1_epi16(16);
    __m128i* result = (__m128i*)(out.data());

    if (border)
      {
      i4 = i2;
      i3 = i1;
      i2 = i0;
      while (result != result_v)
        {
        __m128i ilo, ihi;
        unpack_8bit_to_16bit(*i2, ihi, ilo);
        ihi = _mm_mullo_epi16(ihi, sixes);
        ilo = _mm_mullo_epi16(ilo, sixes);
        *result = _mm_add_epi16(*result, ihi);
        *(result + 1) = _mm_add_epi16(*(result + 1), ilo);
        unpack_8bit_to_16bit(*i3, ihi, ilo);
        ihi = _mm_mullo_epi16(ihi, fours);
        ilo = _mm_mullo_epi16(ilo, fours);
        *result = _mm_add_epi16(*result, ihi);
        *(result + 1) = _mm_add_epi16(*(result + 1), ilo);
        unpack_8bit_to_16bit(*i4, ihi, ilo);
        *result = _mm_add_epi16(*result, ihi);
        *(result + 1) = _mm_add_epi16(*(result + 1), ilo);
        result += 2;
        ++i2; ++i3; ++i4;
        }
      }
    else
      {
      __m128i* i5 = (__m128i*)(im.data());
      while (result != result_v)
        {
        __m128i ilo, ihi;
        unpack_8bit_to_16bit(*i5, ihi, ilo);
        *result = _mm_mullo_epi16(ihi, sixteen);
        ++result;
        *result = _mm_mullo_epi16(ilo, sixteen);
        ++result;
        ++i5;
        }
      }

    for (; i4 != end_input; i0++, i1++, i2++, i3++, i4++, result_v += 2)
      {
      __m128i ilo, ihi;
      unpack_8bit_to_16bit(*i0, ihi, ilo);
      *result_v = _mm_add_epi16(*result_v, ihi);
      *(result_v + 1) = _mm_add_epi16(*(result_v + 1), ilo);
      unpack_8bit_to_16bit(*i1, ihi, ilo);
      ihi = _mm_mullo_epi16(ihi, fours);
      ilo = _mm_mullo_epi16(ilo, fours);
      *result_v = _mm_add_epi16(*result_v, ihi);
      *(result_v + 1) = _mm_add_epi16(*(result_v + 1), ilo);
      unpack_8bit_to_16bit(*i2, ihi, ilo);
      ihi = _mm_mullo_epi16(ihi, sixes);
      ilo = _mm_mullo_epi16(ilo, sixes);
      *result_v = _mm_add_epi16(*result_v, ihi);
      *(result_v + 1) = _mm_add_epi16(*(result_v + 1), ilo);
      unpack_8bit_to_16bit(*i3, ihi, ilo);
      ihi = _mm_mullo_epi16(ihi, fours);
      ilo = _mm_mullo_epi16(ilo, fours);
      *result_v = _mm_add_epi16(*result_v, ihi);
      *(result_v + 1) = _mm_add_epi16(*(result_v + 1), ilo);
      unpack_8bit_to_16bit(*i4, ihi, ilo);
      *result_v = _mm_add_epi16(*result_v, ihi);
      *(result_v + 1) = _mm_add_epi16(*(result_v + 1), ilo);
      }
    if (border)
      {
      while (i2 != end_input)
        {
        __m128i ilo, ihi;
        unpack_8bit_to_16bit(*i2, ihi, ilo);
        ihi = _mm_mullo_epi16(ihi, sixes);
        ilo = _mm_mullo_epi16(ilo, sixes);
        *result_v = _mm_add_epi16(*result_v, ihi);
        *(result_v + 1) = _mm_add_epi16(*(result_v + 1), ilo);
        unpack_8bit_to_16bit(*i1, ihi, ilo);
        ihi = _mm_mullo_epi16(ihi, fours);
        ilo = _mm_mullo_epi16(ilo, fours);
        *result_v = _mm_add_epi16(*result_v, ihi);
        *(result_v + 1) = _mm_add_epi16(*(result_v + 1), ilo);
        unpack_8bit_to_16bit(*i0, ihi, ilo);
        *result_v = _mm_add_epi16(*result_v, ihi);
        *(result_v + 1) = _mm_add_epi16(*(result_v + 1), ilo);
        result_v += 2;
        ++i2; ++i1; ++i0;
        }
      }
    else
      {
      while (i2 != end_input)
        {
        __m128i ilo, ihi;
        unpack_8bit_to_16bit(*i2, ihi, ilo);
        *result_v = _mm_mullo_epi16(ihi, sixteen);
        ++result_v;
        *result_v = _mm_mullo_epi16(ilo, sixteen);
        ++result_v;
        ++i2;
        }
      }
    }


  void convolve_row_14641_div_256_sse(image<uint8_t>& out, const image<int16_t>& im)
    {
    out = image<uint8_t>(im.width(), im.height());
    const __m128i* i0 = (const __m128i*)(im.data());
    const int16_t* i1 = im.data() + 1;
    const int16_t* i2 = im.data();
    const int16_t* i3 = im.data() + 3;
    const int16_t* i4 = im.data() + 4;
    uint8_t* result = out.data();
    *result = (uint8_t)(*i2 >> 4);
    ++i2; ++result;
    *result = (uint8_t)(*i2 >> 4);
    ++i2; ++result;
    const int16_t* const end_input = im.data() + im.stride() * im.height();
    __m128i result_register_lo = _mm_setzero_si128();
    __m128i result_register_hi = _mm_setzero_si128();
    for (; i4 < end_input; i0 += 1, i1 += 8, i2 += 8, i3 += 8, i4 += 8, result += 16)
      {
      for (int i = 0; i < 2; ++i)
        {
        __m128i* result_register;
        if (i == 0)
          result_register = &result_register_lo;
        else
          result_register = &result_register_hi;
        __m128i i0_register = *i0;
        __m128i i1_register = _mm_loadu_si128((__m128i*)(i1));
        __m128i i2_register = _mm_loadu_si128((__m128i*)(i2));
        __m128i i3_register = _mm_loadu_si128((__m128i*)(i3));
        __m128i i4_register = _mm_loadu_si128((__m128i*)(i4));
        *result_register = _mm_setzero_si128();
        *result_register = _mm_add_epi16(i0_register, *result_register);
        i1_register = _mm_add_epi16(i1_register, i1_register);
        i1_register = _mm_add_epi16(i1_register, i1_register);
        *result_register = _mm_add_epi16(i1_register, *result_register);
        i2_register = _mm_add_epi16(i2_register, i2_register);
        *result_register = _mm_add_epi16(i2_register, *result_register);
        i2_register = _mm_add_epi16(i2_register, i2_register);
        *result_register = _mm_add_epi16(i2_register, *result_register);
        i3_register = _mm_add_epi16(i3_register, i3_register);
        i3_register = _mm_add_epi16(i3_register, i3_register);
        *result_register = _mm_add_epi16(i3_register, *result_register);
        *result_register = _mm_add_epi16(i4_register, *result_register);
        *result_register = _mm_srai_epi16(*result_register, 8);
        if (i == 0)
          {
          i0 += 1;
          i1 += 8;
          i2 += 8;
          i3 += 8;
          i4 += 8;
          }
        }
      result_register_lo = _mm_packs_epi16(result_register_lo, result_register_hi);
      _mm_storeu_si128(((__m128i*)(result)), result_register_lo);
      }
    i2 = im.end() - 2;
    result = out.end() - 2;
    *result = (uint8_t)(*i2 >> 4);
    ++i2; ++result;
    *result = (uint8_t)(*i2 >> 4);
    }


  void gauss_sse(image<uint8_t>& out, const image<uint8_t>& im)
    {
    image<int16_t> temp;
    convolve_col_14641_sse(temp, im);
    convolve_row_14641_div_256_sse(out, temp);
    }

  void split(image<uint8_t>& red, image<uint8_t>& green, image<uint8_t>& blue, const image<uint32_t>& im)
    {
    const uint32_t w = im.width();
    const uint32_t h = im.height();
    red = image<uint8_t>(w, h);
    green = image<uint8_t>(w, h);
    blue = image<uint8_t>(w, h);
    for (uint32_t y = 0; y < h; ++y)
      {
      const uint32_t* p_im = im.data() + y * im.stride();
      uint8_t* p_r = red.data() + y * red.stride();
      uint8_t* p_g = green.data() + y * green.stride();
      uint8_t* p_b = blue.data() + y * blue.stride();
      for (uint32_t x = 0; x < w; ++x, ++p_im, ++p_r, ++p_g, ++p_b)
        {
        uint32_t clr = *p_im;
        *p_r = clr & 255;
        *p_g = (clr >> 8) & 255;
        *p_b = (clr >> 16) & 255;
        }
      }
    }

  image<uint8_t> remove_even_col_even_row_sse(const image<uint8_t>& im)
    {
    image<uint8_t> out(im.width() >> 1, (im.height() + im.height() % 2) >> 1);
    const int w_chunk = im.stride() / 32;
    const uint8_t* i0 = im.data();
    const uint8_t* i1 = im.data() + 16;
    __m128i* result = (__m128i*)(out.data());
    const uint8_t* const end_input = im.data() + im.stride() * im.height();
    int cnt = 1;
    for (; i1 < end_input; i0 += 32, i1 += 32, ++result, ++cnt)
      {
      __m128i i0_register = *(const __m128i*)i0;
      __m128i i1_register = *(const __m128i*)i1;
      __m128i hi, lo;
      split_hi_lo(i0_register, i1_register, hi, lo);
      split_hi_lo(lo, hi, hi, lo);
      split_hi_lo(lo, hi, hi, lo);
      split_hi_lo(lo, hi, hi, lo);
      *result = lo;
      if (cnt == w_chunk)
        {
        cnt = 0;
        i0 += im.width();
        i1 += im.width();
        }
      }
    return out;
    }

  image<uint8_t> add_even_col_even_row_sse(const image<uint8_t>& im)
    {
    image<uint8_t> out(im.width() << 1, im.height() << 1);
    const int w_chunk = im.stride() / 16;
    const uint8_t* i0 = im.data();
    __m128i* result = (__m128i*)(out.data());
    const uint8_t* const end_input = im.data() + im.stride() * im.height();
    int cnt = 1;
    while (i0 != end_input)
      {
      __m128i zero = _mm_setzero_si128();
      __m128i b0 = _mm_unpacklo_epi8(*((const __m128i*)i0), zero);
      __m128i b1 = _mm_unpackhi_epi8(*((const __m128i*)i0), zero);
      *result = b0;
      ++result;
      *result = b1;
      ++result;
      i0 += 16;
      if (cnt == w_chunk)
        {
        cnt = 1;
        result += 2 * w_chunk;
        }
      else
        ++cnt;
      }
    return out;
    }
  }

jtk::image<uint32_t> gauss(const jtk::image<uint32_t>& im)
  {
  jtk::image<uint8_t> red, green, blue;
  split(red, green, blue, im);
  gauss_sse(red, red);
  gauss_sse(green, green);
  gauss_sse(blue, blue);
  return three_gray_to_uint32_t(red, green, blue);
  }

jtk::image<uint32_t> pyramid_down(const jtk::image<uint32_t>& im)
  {
  jtk::image<uint8_t> red, green, blue;
  split(red, green, blue, im);
  gauss_sse(red, red);
  gauss_sse(green, green);
  gauss_sse(blue, blue);
  red = remove_even_col_even_row_sse(red);
  green = remove_even_col_even_row_sse(green);
  blue = remove_even_col_even_row_sse(blue);
  return three_gray_to_uint32_t(red, green, blue);
  }

jtk::image<uint32_t> pyramid_up(const jtk::image<uint32_t>& im)
  {
  jtk::image<uint8_t> red, green, blue;
  split(red, green, blue, im);
  jtk::image<int16_t> temp;
  convolve_col_14641_sse(temp, add_even_col_even_row_sse(red), true);
  for (auto& v : temp)
    v *= 4;
  convolve_row_14641_div_256_sse(red, temp);
  convolve_col_14641_sse(temp, add_even_col_even_row_sse(green), true);
  for (auto& v : temp)
    v *= 4;
  convolve_row_14641_div_256_sse(green, temp);
  convolve_col_14641_sse(temp, add_even_col_even_row_sse(blue), true);
  for (auto& v : temp)
    v *= 4;
  convolve_row_14641_div_256_sse(blue, temp);
  return three_gray_to_uint32_t(red, green, blue);
  }
