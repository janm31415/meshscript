#include "mm.h"
#include <jtk/file_utils.h>
#include <jtk/qbvh.h>

bool read_from_file(mm& morph, const std::string& filename)
  {
  std::string ext = jtk::get_extension(filename);
  if (ext.empty())
    return false;
  std::transform(ext.begin(), ext.end(), ext.begin(), [](char ch) {return (char)::tolower(ch); });
  if (ext == "ssm")
    {
    if (!jtk::read_morphable_model_binary(morph.m, filename.c_str()))
      return false;
    }
  else
    return false;
  morph.coefficients.resize(morph.m.S.rows(), 0.f);
  morph.vertices = jtk::get_vertices(morph.m, morph.coefficients);  
  morph.cs = jtk::get_identity();
  morph.visible = true;
  return true;
  }