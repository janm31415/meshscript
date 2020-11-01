#include "im.h"

#include <stb_image.h>
#include <stb_image_write.h>

#include <jtk/file_utils.h>

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