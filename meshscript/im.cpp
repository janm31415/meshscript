#include "im.h"

#include <stb_image.h>
#include <stb_image_write.h>

#include <jtk/file_utils.h>

#include <algorithm>
#include <iostream>

bool read_from_file(im& i, const std::string& filename)
  {
  int w, h, nr_of_channels;
  unsigned char* im = stbi_load(filename.c_str(), &w, &h, &nr_of_channels, 0);
  if (im)
    {
    i.texture = jtk::image<uint32_t>(w, h);
    for (uint32_t y = 0; y < (uint32_t)h; ++y)
      {
      uint32_t* p_out = i.texture.data() + y * w;
      const unsigned char* p_im = im + y * w * nr_of_channels;
      for (uint32_t x = 0; x < (uint32_t)w; ++x, ++p_out, p_im += nr_of_channels)
        {
        uint32_t color = 0;
        switch (nr_of_channels)
          {
          case 1:
          {
          uint32_t g = *p_im;
          color = 0xff000000 | (g << 16) | (g << 8) | g;
          break;
          }
          case 3:
          {
          uint32_t r = p_im[0];
          uint32_t g = p_im[1];
          uint32_t b = p_im[2];
          color = 0xff000000 | (b << 16) | (g << 8) | r;
          break;
          }
          case 4:
          {
          color = *((const uint32_t*)p_im);
          break;
          }
          }
        *p_out = color;
        }
      }
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
