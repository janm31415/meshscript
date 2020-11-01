#pragma once

#include <jtk/image.h>
#include <string>

struct im
  {
  jtk::image<uint32_t> texture;
  };


bool read_from_file(im& i, const std::string& filename);
bool write_to_file(const im& i, const std::string& filename);

void info(const im& i);