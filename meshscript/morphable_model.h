#pragma once

#include <jtk/mat.h>
#include <jtk/vec.h>

#include <vector>

namespace jtk
  {

  struct morphable_model
    {
    matf average;
    matf U;
    matf V;
    matf S;

    std::vector<vec3<uint32_t>> triangles;
    };

  void swap(morphable_model& left, morphable_model& right);
  std::vector<vec3<float>> get_vertices(const morphable_model& mm, const std::vector<float>& coefficients);
  std::vector<float> get_basic_shape(const morphable_model& mm, uint32_t index);
  bool read_morphable_model_binary(morphable_model& mm, const char* filename);
  bool write_morphable_model_binary(const morphable_model& mm, const char* filename);
  float sigma(const morphable_model& mm, uint64_t index);

  inline void swap(morphable_model& left, morphable_model& right)
    {
    left.average.swap(right.average);
    left.U.swap(right.U);
    left.V.swap(right.V);
    left.S.swap(right.S);
    std::swap(left.triangles, right.triangles);
    }

  inline std::vector<vec3<float>> get_vertices(const morphable_model& mm, const std::vector<float>& coefficients)
    {
    uint64_t dim1 = mm.U.rows() / 3;
    uint64_t dim2 = (uint64_t)coefficients.size();
    if (dim2 > mm.U.cols())
      dim2 = mm.U.cols();
    std::vector<vec3<float>> vertices(dim1);
    for (uint64_t i = 0; i < dim1; ++i)
      {
      vec3<float> pt(mm.average(i * 3), mm.average(i * 3 + 1), mm.average(i * 3 + 2));
      for (uint64_t j = 0; j < dim2; ++j)
        {
        pt[0] += coefficients[j] * mm.U(i * 3, j);
        pt[1] += coefficients[j] * mm.U(i * 3 + 1, j);
        pt[2] += coefficients[j] * mm.U(i * 3 + 2, j);
        }
      vertices[i] = pt;
      }
    return vertices;
    }

  inline std::vector<float> get_basic_shape(const morphable_model& mm, uint32_t index)
    {
    std::vector<float> coeff(mm.S.rows(), 0.f);
    for (size_t i = 0; i < coeff.size(); ++i)
      {
      coeff[i] = mm.S(i) * mm.V(index, i);
      }
    return coeff;
    }


  inline bool read_morphable_model_binary(morphable_model& mm, const char* filename)
    {
    FILE* inputfile;

    inputfile = fopen(filename, "rb");

    if (!inputfile)
      return false;

    // read average vector
    uint32_t sz;
    fread(&sz, sizeof(uint32_t), 1, inputfile);
    mm.average.resize(sz, 1);
    fread((void*)mm.average.data(), sizeof(float), sz, inputfile);

    // read U matrix (stored row-major)
    uint32_t d1, d2;
    fread(&d1, sizeof(uint32_t), 1, inputfile);
    fread(&d2, sizeof(uint32_t), 1, inputfile);
    mm.U.resize(d1, d2);
    fread((void*)mm.U.data(), sizeof(float), d1*d2, inputfile);

    // read V matrix (stored row-major)
    fread(&d1, sizeof(uint32_t), 1, inputfile);
    fread(&d2, sizeof(uint32_t), 1, inputfile);
    mm.V.resize(d1, d2);
    fread((void*)mm.V.data(), sizeof(float), d1*d2, inputfile);

    // read S vector
    fread(&sz, sizeof(uint32_t), 1, inputfile);
    mm.S.resize(sz, 1);
    fread((void*)mm.S.data(), sizeof(float), sz, inputfile);

    // read triangles vector
    fread(&sz, sizeof(uint32_t), 1, inputfile);
    assert(sz % 3 == 0);
    mm.triangles.resize(sz / 3);
    fread((void*)mm.triangles.data(), sizeof(uint32_t), sz, inputfile);

    fclose(inputfile);
    return true;
    }

  inline bool write_morphable_model_binary(const morphable_model& mm, const char* filename)
    {
    FILE* outputfile;
    outputfile = fopen(filename, "wb");

    if (!outputfile)
      return false;

    // write average vector
    uint32_t sz = (uint32_t)mm.average.rows();
    fwrite(&sz, sizeof(uint32_t), 1, outputfile);
    fwrite((void*)mm.average.data(), sizeof(float), sz, outputfile);

    // write U matrix (row-major)
    uint32_t d1 = (uint32_t)mm.U.rows();
    uint32_t d2 = (uint32_t)mm.U.cols();
    fwrite(&d1, sizeof(uint32_t), 1, outputfile);
    fwrite(&d2, sizeof(uint32_t), 1, outputfile);   
    fwrite((void*)mm.U.data(), sizeof(float), d1*d2, outputfile);

    // write V matrix (row-major)
    d1 = (uint32_t)mm.V.rows();
    d2 = (uint32_t)mm.V.cols();
    fwrite(&d1, sizeof(uint32_t), 1, outputfile);
    fwrite(&d2, sizeof(uint32_t), 1, outputfile);
    fwrite((void*)mm.V.data(), sizeof(float), d1*d2, outputfile);

    // write S vector
    sz = (uint32_t)mm.S.rows();
    fwrite(&sz, sizeof(uint32_t), 1, outputfile);
    fwrite((void*)mm.S.data(), sizeof(float), sz, outputfile);

    // read triangles vector
    sz = (uint32_t)mm.triangles.size() * 3;
    fwrite(&sz, sizeof(uint32_t), 1, outputfile);
    fwrite((void*)mm.triangles.data(), sizeof(uint32_t), sz, outputfile);

    fclose(outputfile);
    return true;
    }

  inline float sigma(const morphable_model& mm, uint64_t index)
    {
    float denom = float(mm.S.rows());
    if (denom > 1)
      denom -= 1.f;
    return mm.S(index) / std::sqrt(denom);
    }

  } // namespace jtk
