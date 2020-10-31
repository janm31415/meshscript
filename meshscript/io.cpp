#include "io.h"
#include <string.h>

#include "rply.h"

#include <trico/trico/trico.h>

#include <iostream>
#include <jtk/file_utils.h>

#include <stb_image.h>
#include <stb_image_write.h>

static int read_vec3_coord(p_ply_argument argument)
  {
  float** p_vertices;
  ply_get_argument_user_data(argument, (void**)&p_vertices, NULL);
  double val = ply_get_argument_value(argument);
  *(*p_vertices) = (float)val;
  (*p_vertices) += 3;
  return 1;
  }

static int read_color(p_ply_argument argument)
  {
  uint8_t** p_color;
  ply_get_argument_user_data(argument, (void**)&p_color, NULL);
  double val = ply_get_argument_value(argument);
  *(*p_color) = (uint8_t)val;
  (*p_color) += 4;
  return 1;
  }

static int read_face(p_ply_argument argument)
  {
  uint32_t** p_triangle;
  ply_get_argument_user_data(argument, (void**)&p_triangle, NULL);
  long length, value_index;
  ply_get_argument_property(argument, NULL, &length, &value_index);
  if (value_index >= 0 && value_index < 3)
    {
    double val = ply_get_argument_value(argument);
    *(*p_triangle) = (uint32_t)val;
    ++(*p_triangle);
    }
  return 1;
  }

static int read_texcoord(p_ply_argument argument)
  {
  float** p_uv;
  ply_get_argument_user_data(argument, (void**)&p_uv, NULL);
  long length, value_index;
  ply_get_argument_property(argument, NULL, &length, &value_index);
  if (value_index >= 0 && value_index < 6)
    {
    double val = ply_get_argument_value(argument);
    *(*p_uv) = (float)val;
    ++(*p_uv);
    }
  if (value_index == (length - 1) && (length != 6))
    {
    for (long j = length; j < 6; ++j)
      {
      *(*p_uv) = (float)0.f;
      ++(*p_uv);
      }
    }
  return 1;
  }

bool read_ply(const char* filename, std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<jtk::vec2<float>>>& uv)
  {
  vertices.clear();
  normals.clear();
  clrs.clear();
  triangles.clear();
  uv.clear();

  p_ply ply = ply_open(filename, NULL, 0, NULL);
  if (!ply) 
    return false;
  if (!ply_read_header(ply)) 
    return false;

  float* p_vertex_pointer_x = NULL;
  float* p_vertex_pointer_y = NULL;
  float* p_vertex_pointer_z = NULL;

  long nvertices_x = ply_set_read_cb(ply, "vertex", "x", read_vec3_coord, (void*)(&p_vertex_pointer_x), 0);
  long nvertices_y = ply_set_read_cb(ply, "vertex", "y", read_vec3_coord, (void*)(&p_vertex_pointer_y), 0);
  long nvertices_z = ply_set_read_cb(ply, "vertex", "z", read_vec3_coord, (void*)(&p_vertex_pointer_z), 0);

  if (nvertices_x != nvertices_y)
    return false;
  if (nvertices_x != nvertices_z)
    return false;

  if (nvertices_x > 0)
    vertices.resize(nvertices_x);
  p_vertex_pointer_x = (float*)vertices.data();
  p_vertex_pointer_y = p_vertex_pointer_x + 1;
  p_vertex_pointer_z = p_vertex_pointer_x + 2;

  float* p_normal_pointer_x = NULL;
  float* p_normal_pointer_y = NULL;
  float* p_normal_pointer_z = NULL;

  long nnormals_x = ply_set_read_cb(ply, "vertex", "nx", read_vec3_coord, (void*)(&p_normal_pointer_x), 0);
  long nnormals_y = ply_set_read_cb(ply, "vertex", "ny", read_vec3_coord, (void*)(&p_normal_pointer_y), 0);
  long nnormals_z = ply_set_read_cb(ply, "vertex", "nz", read_vec3_coord, (void*)(&p_normal_pointer_z), 0);

  if (nnormals_x != nnormals_y)
    return false;
  if (nnormals_x != nnormals_z)
    return false;

  if (nnormals_x > 0)
    normals.resize(nnormals_x);
  p_normal_pointer_x = (float*)normals.data();
  p_normal_pointer_y = p_normal_pointer_x + 1;
  p_normal_pointer_z = p_normal_pointer_x + 2;

  uint8_t* p_red = NULL;
  uint8_t* p_green = NULL;
  uint8_t* p_blue = NULL;
  uint8_t* p_alpha = NULL;

  long nred = ply_set_read_cb(ply, "vertex", "red", read_color, (void*)(&p_red), 0);
  long ngreen = ply_set_read_cb(ply, "vertex", "green", read_color, (void*)(&p_green), 0);
  long nblue = ply_set_read_cb(ply, "vertex", "blue", read_color, (void*)(&p_blue), 0);
  long nalpha = ply_set_read_cb(ply, "vertex", "alpha", read_color, (void*)(&p_alpha), 0);

  if (nred == 0)
    nred = ply_set_read_cb(ply, "vertex", "r", read_color, (void*)(&p_red), 0);
  if (ngreen == 0)
    ngreen = ply_set_read_cb(ply, "vertex", "g", read_color, (void*)(&p_green), 0);
  if (nblue == 0)
    nblue = ply_set_read_cb(ply, "vertex", "b", read_color, (void*)(&p_blue), 0);
  if (nalpha == 0)
    nalpha = ply_set_read_cb(ply, "vertex", "a", read_color, (void*)(&p_alpha), 0);

  if (nred == 0)    
    nred = ply_set_read_cb(ply, "vertex", "diffuse_red", read_color, (void*)(&p_red), 0);
  if (ngreen == 0)
    ngreen = ply_set_read_cb(ply, "vertex", "diffuse_green", read_color, (void*)(&p_green), 0);
  if (nblue == 0)
    nblue = ply_set_read_cb(ply, "vertex", "diffuse_blue", read_color, (void*)(&p_blue), 0);
  if (nalpha == 0)
    nalpha = ply_set_read_cb(ply, "vertex", "diffuse_alpha", read_color, (void*)(&p_alpha), 0);      

  if (nred > 0 || ngreen > 0 || nblue > 0 || nalpha > 0)
    clrs.resize(nred, 0xffffffff);

  p_red = (uint8_t*)clrs.data();
  p_green = p_red + 1;
  p_blue = p_red + 2;
  p_alpha = p_red + 3;

  uint32_t* p_tria_index = NULL;

  long ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", read_face, (void*)(&p_tria_index), 0);
  if (ntriangles == 0)
    ntriangles = ply_set_read_cb(ply, "face", "vertex_index", read_face, (void*)(&p_tria_index), 0);

  if (ntriangles > 0)
    triangles.resize(ntriangles);

  p_tria_index = (uint32_t*)triangles.data();

  float* p_uv = NULL;

  long ntexcoords = ply_set_read_cb(ply, "face", "texcoord", read_texcoord, (void*)(&p_uv), 0);

  if (ntexcoords > 0)
    uv.resize(ntexcoords);

  p_uv = (float*)uv.data();

  if (!ply_read(ply))
    return false;

  ply_close(ply);

  return true;
  }

/*
extern "C"
  {
#include "ply.h"
  }

namespace
  {

  struct ply_vert
    {
    float x, y, z, nx, ny, nz;
    unsigned char r, g, b, a;
    };

  struct ply_tria
    {
    unsigned char nverts;
    int* verts;
    unsigned char ntexcoords;
    float* texcoords;
    };

  PlyProperty ply_vert_props[] = { 
  {"x", PLY_FLOAT, PLY_FLOAT, offsetof(ply_vert, x), 0, 0, 0, 0},
  {"y", PLY_FLOAT, PLY_FLOAT, offsetof(ply_vert, y), 0, 0, 0, 0},
  {"z", PLY_FLOAT, PLY_FLOAT, offsetof(ply_vert, z), 0, 0, 0, 0},
  {"nx", PLY_FLOAT, PLY_FLOAT, offsetof(ply_vert, nx), 0, 0, 0, 0},
  {"ny", PLY_FLOAT, PLY_FLOAT, offsetof(ply_vert, ny), 0, 0, 0, 0},
  {"nz", PLY_FLOAT, PLY_FLOAT, offsetof(ply_vert, nz), 0, 0, 0, 0},
  {"diffuse_red", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, r), 0, 0, 0, 0},
  {"diffuse_green", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, g), 0, 0, 0, 0},
  {"diffuse_blue", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, b), 0, 0, 0, 0},
  {"diffuse_alpha", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, b), 0, 0, 0, 0},
  {"red", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, r), 0, 0, 0, 0},
  {"green", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, g), 0, 0, 0, 0},
  {"blue", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, b), 0, 0, 0, 0},
  {"alpha", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, a), 0, 0, 0, 0},
  {"r", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, r), 0, 0, 0, 0},
  {"g", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, g), 0, 0, 0, 0},
  {"b", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, b), 0, 0, 0, 0},
  {"a", PLY_UCHAR, PLY_UCHAR, offsetof(ply_vert, a), 0, 0, 0, 0}
    };

  PlyProperty ply_face_props[] = {
      {"vertex_indices", PLY_INT, PLY_INT, offsetof(ply_tria, verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(ply_tria, nverts)},
      {"vertex_index", PLY_INT, PLY_INT, offsetof(ply_tria, verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(ply_tria, nverts)},
      {"texcoord", PLY_FLOAT, PLY_FLOAT, offsetof(ply_tria, texcoords), 1, PLY_UCHAR, PLY_UCHAR, offsetof(ply_tria, ntexcoords)}
    };

  }

bool read_ply(const char* filename, std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<jtk::vec2<float>>>& uv)
  {
  vertices.clear();
  normals.clear();
  clrs.clear();
  triangles.clear();
  uv.clear();

  bool has_vertices = false;
  bool has_color = false;
  bool has_faces = false;
  bool has_texcoords = false;
  bool has_normals = false;

  PlyFile* ply;
  int nr_elems;
  char **elist;
  int file_type;
  float version;

  ply = ply_open_for_reading(filename, &nr_elems, &elist, &file_type, &version);
  if (!ply)
    return false;


  char* elem_name;
  int num_elems;
  int i, j;
  PlyProperty** plist;
  int nr_props;

  bool result = true;

  for (i = 0; i < nr_elems; ++i)
    {
    elem_name = elist[i];
    plist = ply_get_element_description(ply, elem_name, &num_elems, &nr_props);
    if (!plist)
      {
      for (i = 0; i < nr_elems; i++) {
        free(ply->elems[i]->name);
        free(ply->elems[i]->store_prop);
        for (j = 0; j < ply->elems[i]->nprops; j++) {
          free(ply->elems[i]->props[j]->name);
          free(ply->elems[i]->props[j]);
          }
        free(ply->elems[i]->props);
        }
      for (i = 0; i < nr_elems; i++) { free(ply->elems[i]); }
      free(ply->elems);
      for (i = 0; i < ply->num_comments; i++) { free(ply->comments[i]); }
      free(ply->comments);
      for (i = 0; i < ply->num_obj_info; i++) { free(ply->obj_info[i]); }
      free(ply->obj_info);
      ply_free_other_elements(ply->other_elems);

      for (i = 0; i < nr_elems; i++) { free(elist[i]); }
      free(elist);
      ply_close(ply);
      return false;
      }
    if (equal_strings("vertex", elem_name))
      {
      for (int p = 0; p < nr_props; ++p)
        {
        for (int vp = 0; vp < 18; ++vp)
          {
          if (strcmp(plist[p]->name, ply_vert_props[vp].name) == 0)
            {
            if (vp >= 0 && vp < 3)
              has_vertices = true;
            else if (vp >= 3 && vp < 6)
              has_normals = true;
            else if (vp >= 6 && vp < 18)
              has_color = true;
            ply_get_property(ply, elem_name, &ply_vert_props[vp]);
            if (vp >= 0 && vp < 3)
              {
              if (ply_vert_props[vp].external_type != PLY_FLOAT)
                {
                fprintf(stderr, "Warning: vertices should be of type float\n");
                }
              }
            else if (vp >= 3 && vp < 6)
              {
              if (ply_vert_props[vp].external_type != PLY_FLOAT)
                {
                fprintf(stderr, "Warning: normals should be of type float\n");
                }
              }
            else if (vp >= 6 && vp < 18)
              {
              if (ply_vert_props[vp].external_type != PLY_UCHAR)
                {
                fprintf(stderr, "Warning: colors should be of type uchar\n");
                }
              }
            break;
            }
          }

        }

      if (has_vertices)
        vertices.resize(num_elems);
      if (has_color)
        clrs.resize(num_elems);
      if (has_normals)
        normals.resize(num_elems);

      for (j = 0; j < num_elems; j++)
        {
        ply_vert v;
        v.a = 0xff;
        v.r = v.g = v.b = 0;
        v.x = v.y = v.z = 0.f;
        v.nx = v.ny = v.nz = 0.f;
        ply_get_element(ply, (void *)&v);
        if (has_vertices)
          {
          vertices[j][0] = v.x;
          vertices[j][1] = v.y;
          vertices[j][2] = v.z;
          }
        if (has_normals)
          {
          normals[j][0] = v.nx;
          normals[j][1] = v.ny;
          normals[j][2] = v.nz;
          }
        if (has_color)
          {
          uint32_t r = (uint32_t)v.r;
          uint32_t g = (uint32_t)v.g;
          uint32_t b = (uint32_t)v.b;
          uint32_t a = (uint32_t)v.a;
          clrs[j] = (a << 24) | (b << 16) | (g << 8) | r;
          }
        }
      }
    if (equal_strings("face", elem_name))
      {
      for (int p = 0; p < nr_props; ++p)
        {
        for (int vp = 0; vp < 3; ++vp)
          {
          if (strcmp(plist[p]->name, ply_face_props[vp].name) == 0)
            {
            if (vp == 0 || vp == 1)
              has_faces = true;
            if (vp == 2)
              has_texcoords = true;
            ply_get_property(ply, elem_name, &ply_face_props[vp]);
            if (vp == 0 || vp == 1)
              {
              if (ply_face_props[vp].external_type != PLY_INT)
                {
                fprintf(stderr, "Warning: vertex indices should be of type int\n");
                }
              if (ply_face_props[vp].count_external != PLY_UCHAR)
                {
                fprintf(stderr, "Warning: vertex index list length should be of type uchar\n");
                }
              }
            if (vp == 2)
              {
              if (ply_face_props[vp].external_type != PLY_FLOAT)
                {
                fprintf(stderr, "Warning: texture coordinates should be of type float\n");
                }
              if (ply_face_props[vp].count_external != PLY_UCHAR)
                {
                fprintf(stderr, "Warning: texture coordinates list length should be of type uchar\n");
                }
              }
            break;
            }
          }
        }
      if (has_faces)
        triangles.resize(num_elems);
      if (has_texcoords)
        uv.resize(num_elems);

      for (j = 0; j < num_elems; j++)
        {
        ply_tria t;
        ply_get_element(ply, (void *)&t);

        if (has_faces)
          {
          if (t.nverts == 3)
            {
            triangles[j][0] = t.verts[0];
            triangles[j][1] = t.verts[1];
            triangles[j][2] = t.verts[2];
            }
          else
            {
            triangles[j][0] = 0;
            triangles[j][1] = 0;
            triangles[j][2] = 0;
            fprintf(stderr, "Warning: only triangular faces are supported\n");
            result = false;
            }
          free(t.verts);
          }
        if (has_texcoords)
          {
          if (t.ntexcoords == 6)
            {
            uv[j][0][0] = t.texcoords[0];
            uv[j][0][1] = t.texcoords[1];
            uv[j][1][0] = t.texcoords[2];
            uv[j][1][1] = t.texcoords[3];
            uv[j][2][0] = t.texcoords[4];
            uv[j][2][1] = t.texcoords[5];
            }
          else
            {
            uv[j][0][0] = uv[j][0][1] = 0.f;
            uv[j][1][0] = uv[j][1][1] = 0.f;
            uv[j][2][0] = uv[j][2][1] = 0.f;
            fprintf(stderr, "Warning: 6 texture coordinates per triangular face are expected\n");
            result = false;
            }
          free(t.texcoords);
          }
        }
      }
    for (j = 0; j < nr_props; j++)
      {
      free(plist[j]->name);
      free(plist[j]);
      }
    free(plist);
    }  // for each type of element

  for (i = 0; i < nr_elems; i++)
    {
    free(ply->elems[i]->name);
    free(ply->elems[i]->store_prop);
    for (j = 0; j < ply->elems[i]->nprops; j++) {
      free(ply->elems[i]->props[j]->name);
      free(ply->elems[i]->props[j]);
      }
    if (ply->elems[i]->props && ply->elems[i]->nprops)
      {
      free(ply->elems[i]->props);
      }
    }
  for (i = 0; i < nr_elems; i++) { free(ply->elems[i]); }
  free(ply->elems);
  for (i = 0; i < ply->num_comments; i++) { free(ply->comments[i]); }
  free(ply->comments);
  for (i = 0; i < ply->num_obj_info; i++) { free(ply->obj_info[i]); }
  free(ply->obj_info);
  ply_free_other_elements(ply->other_elems);


  for (i = 0; i < nr_elems; i++) { free(elist[i]); }
  free(elist);
  ply_close(ply);
  return result;
  }

  */


bool write_ply(const char* filename, const std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<float>>& normals, const std::vector<uint32_t>& clrs, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<jtk::vec2<float>>>& uv)
  {
  if (vertices.empty())
    return false;
  FILE* fp = fopen(filename, "wb");

  if (!fp)
    return false;

  fprintf(fp, "ply\n");
  int n = 1;
  if (*(char *)&n == 1)
    fprintf(fp, "format binary_little_endian 1.0\n");
  else
    fprintf(fp, "format binary_big_endian 1.0\n");

  fprintf(fp, "element vertex %d\n", (uint32_t)vertices.size());
  fprintf(fp, "property float x\n");
  fprintf(fp, "property float y\n");
  fprintf(fp, "property float z\n");

  if (!normals.empty())
    {
    fprintf(fp, "property float nx\n");
    fprintf(fp, "property float ny\n");
    fprintf(fp, "property float nz\n");
    }

  if (!clrs.empty())
    {
    fprintf(fp, "property uchar red\n");
    fprintf(fp, "property uchar green\n");
    fprintf(fp, "property uchar blue\n");
    fprintf(fp, "property uchar alpha\n");
    }

  if (!triangles.empty())
    {
    fprintf(fp, "element face %d\n", (uint32_t)triangles.size());
    fprintf(fp, "property list uchar int vertex_indices\n");
    if (!uv.empty())
      fprintf(fp, "property list uchar float texcoord\n");
    }
  fprintf(fp, "end_header\n");

  for (uint32_t i = 0; i < (uint32_t)vertices.size(); ++i)
    {
    fwrite((float*)vertices.data() + 3 * i, sizeof(float), 3, fp);
    if (!normals.empty())
      fwrite((float*)normals.data() + 3 * i, sizeof(float), 3, fp);
    if (!clrs.empty())
      fwrite((uint32_t*)clrs.data() + i, sizeof(uint32_t), 1, fp);
    }
  const unsigned char tria_size = 3;
  const unsigned char texcoord_size = 6;
  for (uint32_t i = 0; i < (uint32_t)triangles.size(); ++i)
    {
    fwrite(&tria_size, 1, 1, fp);
    fwrite((uint32_t*)triangles.data() + 3 * i, sizeof(uint32_t), 3, fp);
    if (!uv.empty())
      {
      fwrite(&texcoord_size, 1, 1, fp);
      fwrite((float*)uv.data() + 6 * i, sizeof(float), 6, fp);
      }
    }
  fclose(fp);
  return 1;
  }


bool read_trc(const char* filename, std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<jtk::vec2<float>>>& uv)
  {
  long long size = jtk::file_size(filename);
  if (size < 0)
    {
    std::cout << "There was an error reading file " << filename << std::endl;
    return false;
    }
  FILE* f = fopen(filename, "rb");
  if (!f)
    {
    std::cout << "Cannot open file: " << filename << std::endl;
    return false;
    }
  char* buffer = (char*)malloc(size);
  long long fl = (long long)fread(buffer, 1, size, f);
  if (fl != size)
    {
    std::cout << "There was an error reading file " << filename << std::endl;
    fclose(f);
    return false;
    }
  fclose(f);

  void* arch = trico_open_archive_for_reading((const uint8_t*)buffer, size);
  if (!arch)
    {
    std::cout << "The input file " << filename << " is not a trico archive." << std::endl;
    return false;
    }

  enum trico_stream_type st = trico_get_next_stream_type(arch);
  while (!st == trico_empty)
    {
    switch (st)
      {
      case trico_vertex_float_stream:
      {
      vertices.resize(trico_get_number_of_vertices(arch));
      float* vert = (float*)vertices.data();
      if (!trico_read_vertices(arch, &vert))
        {
        std::cout << "Something went wrong reading the vertices" << std::endl;
        trico_close_archive(arch);
        free(buffer);
        return false;
        }
      break;
      }
      case trico_triangle_uint32_stream:
      {
      triangles.resize(trico_get_number_of_triangles(arch));
      uint32_t* tria = (uint32_t*)triangles.data();
      if (!trico_read_triangles(arch, &tria))
        {
        std::cout << "Something went wrong reading the triangles" << std::endl;
        trico_close_archive(arch);
        free(buffer);
        return false;
        }
      break;
      }
      case trico_vertex_color_stream:
      {
      clrs.resize(trico_get_number_of_colors(arch));
      uint32_t* vertex_colors = (uint32_t*)clrs.data();
      if (!trico_read_vertex_colors(arch, &vertex_colors))
        {
        std::cout << "Something went wrong reading the vertex colors" << std::endl;
        trico_close_archive(arch);
        free(buffer);
        return false;
        }
      break;      
      }
      case trico_uv_per_triangle_float_stream:
      {
      uv.resize(trico_get_number_of_uvs(arch));
      float* uvs = (float*)uv.data();
      if (!trico_read_uv_per_triangle(arch, &uvs))
        {
        std::cout << "Something went wrong reading the uv coordinates" << std::endl;
        trico_close_archive(arch);
        free(buffer);
        return false;
        }
      break;
      }
      case trico_vertex_normal_float_stream:
      {
      normals.resize(trico_get_number_of_normals(arch));
      float* norm = (float*)normals.data();
      if (!trico_read_vertex_normals(arch, &norm))
        {
        std::cout << "Something went wrong reading the normals" << std::endl;
        trico_close_archive(arch);
        free(buffer);
        return false;
        }
      break;
      }
      default:
      {
      trico_skip_next_stream(arch);
      break;
      }
      }

    st = trico_get_next_stream_type(arch);
    }


  trico_close_archive(arch);
  free(buffer);

  return true;
  }


bool write_trc(const char* filename, const std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<float>>& normals, const std::vector<uint32_t>& clrs, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<jtk::vec2<float>>>& uv)
  {
  if (vertices.empty())
    return false;
  void* arch = trico_open_archive_for_writing(1024 * 1024);
  if (!trico_write_vertices(arch, (float*)vertices.data(), (uint32_t)vertices.size()))
    {
    std::cout << "Something went wrong when writing the vertices\n";
    return false;
    }
  if (!clrs.empty() && !trico_write_vertex_colors(arch, (uint32_t*)clrs.data(), (uint32_t)clrs.size()))
    {
    std::cout << "Something went wrong when writing the vertex colors\n";
    return false;
    }
  if (!normals.empty() && !trico_write_vertex_normals(arch, (float*)normals.data(), (uint32_t)normals.size()))
    {
    std::cout << "Something went wrong when writing the normals\n";
    return false;
    }
  if (!triangles.empty() && !trico_write_triangles(arch, (uint32_t*)triangles.data(), (uint32_t)triangles.size()))
    {
    std::cout << "Something went wrong when writing the triangles\n";
    return false;
    }
  if (!uv.empty() && !trico_write_uv_per_triangle(arch, (float*)uv.data(), (uint32_t)uv.size()))
    {
    std::cout << "Something went wrong when writing the uv positions\n";
    return false;
    }

  FILE* f = fopen(filename, "wb");
  if (!f)
    {
    std::cout << "Cannot write to file " << filename << std::endl;
    return false;
    }

  fwrite((const void*)trico_get_buffer_pointer(arch), trico_get_size(arch), 1, f);
  fclose(f);

  trico_close_archive(arch);

  return true;
  }

namespace
  {

  bool read_texture_filename_from_mtl(std::string& texture_file, const char* filename)
    {
    FILE* f = nullptr;
    f = fopen(filename, "r");
    if (!f)
      return false;
    char buffer[256];
    while (fgets(buffer, 256, f) != nullptr)
      {
      int first_non_whitespace_index = 0;
      while (first_non_whitespace_index < 250 && (buffer[first_non_whitespace_index] == ' ' || buffer[first_non_whitespace_index] == '\t'))
        ++first_non_whitespace_index;
      if (buffer[first_non_whitespace_index + 0] == 'm' && buffer[first_non_whitespace_index + 1] == 'a' && buffer[first_non_whitespace_index + 2] == 'p' && buffer[first_non_whitespace_index + 3] == '_' && buffer[first_non_whitespace_index + 4] == 'K' && buffer[first_non_whitespace_index + 5] == 'd')
        {
        char texture[256];
        auto scan_err = sscanf(buffer + first_non_whitespace_index, "map_Kd %s\n", texture);
        scan_err;
        texture_file = std::string(texture);
        fclose(f);
        return true;
        }
      }
    fclose(f);
    return false;
    }

  }

bool read_obj(const char* filename, std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<float>>& normals, std::vector<uint32_t>& clrs, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<jtk::vec3<jtk::vec2<float>>>& uv, jtk::image<uint32_t>& texture)
  {
  using namespace jtk;
  std::string mtl_filename;
  std::vector<vec3<float>>().swap(vertices);
  std::vector<vec3<float>>().swap(normals);
  std::vector<uint32_t>().swap(clrs);
  std::vector<vec3<uint32_t>>().swap(triangles);
  std::vector<vec3<vec2<float>>>().swap(uv);
  texture = jtk::image<uint32_t>();

  FILE* f = nullptr;
  f = fopen(filename, "r");
  if (!f)
    return false;

  std::vector<vec2<float>> tex;
  std::vector<vec3<uint32_t>> tria_uv;

  char buffer[256];
  while (fgets(buffer, 256, f) != nullptr)
    {
    int first_non_whitespace_index = 0;
    while (first_non_whitespace_index < 250 && (buffer[first_non_whitespace_index] == ' ' || buffer[first_non_whitespace_index] == '\t'))
      ++first_non_whitespace_index;

    if (buffer[first_non_whitespace_index + 0] == 'm' && buffer[first_non_whitespace_index + 1] == 't' && buffer[first_non_whitespace_index + 2] == 'l' && buffer[first_non_whitespace_index + 3] == 'l' && buffer[first_non_whitespace_index + 4] == 'i' && buffer[first_non_whitespace_index + 5] == 'b')
      {
      char mtlfile[256];
      auto scan_err = sscanf(buffer + first_non_whitespace_index, "mtllib %s\n", mtlfile);
      if (scan_err == 1)
        {
        mtl_filename = std::string(mtlfile);
        }
      }
    if (buffer[first_non_whitespace_index + 0] == 'v' && buffer[first_non_whitespace_index + 1] == ' ')
      {
      float x, y, z;
      auto err = sscanf(buffer + first_non_whitespace_index, "v %f %f %f\n", &x, &y, &z);
      if (err != 3)
        {
        fclose(f);
        return false;
        }
      vertices.push_back(vec3<float>(x, y, z));
      }
    else if (buffer[first_non_whitespace_index + 0] == 'v' && buffer[first_non_whitespace_index + 1] == 'n' && buffer[2] == ' ')
      {
      float x, y, z;
      auto err = sscanf(buffer + first_non_whitespace_index, "vn %f %f %f\n", &x, &y, &z);
      if (err != 3)
        {
        fclose(f);
        return false;
        }
      normals.push_back(vec3<float>(x, y, z));
      }
    else if (buffer[first_non_whitespace_index + 0] == 'v' && buffer[first_non_whitespace_index + 1] == 't' && buffer[2] == ' ')
      {
      float x, y;
      auto err = sscanf(buffer + first_non_whitespace_index, "vt %f %f\n", &x, &y);
      if (err != 2)
        {
        fclose(f);
        return false;
        }
      tex.push_back(vec2<float>(x, y));
      }
    else if (buffer[first_non_whitespace_index + 0] == 'f' && buffer[first_non_whitespace_index + 1] == ' ')
      {
      uint32_t t0, t1, t2, t3, v0, v1, v2, v3;
      auto err = sscanf(buffer + first_non_whitespace_index, "f %d/%d %d/%d %d/%d %d/%d\n", &t0, &v0, &t1, &v1, &t2, &v2, &t3, &v3);
      if (err == 6)
        {
        tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
        triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
        }
      else if (err == 8)
        {
        tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
        triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
        tria_uv.push_back(vec3<uint32_t>(v0 - 1, v2 - 1, v3 - 1));
        triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
        }
      else
        {
        err = sscanf(buffer + first_non_whitespace_index, "f %d/%d/ %d/%d/ %d/%d/ %d/%d/\n", &t0, &v0, &t1, &v1, &t2, &v2, &t3, &v3);
        if (err == 6)
          {
          tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
          triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
          }
        else if (err == 8)
          {
          tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
          triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
          tria_uv.push_back(vec3<uint32_t>(v0 - 1, v2 - 1, v3 - 1));
          triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
          }
        else
          {
          err = sscanf(buffer + first_non_whitespace_index, "f %d %d %d %d\n", &t0, &t1, &t2, &t3);
          if (err == 3)
            {
            triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
            }
          else if (err == 4)
            {
            triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
            triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
            }
          else
            {
            uint32_t tx0, tx1, tx2, tx3;
            err = sscanf(buffer + first_non_whitespace_index, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d\n", &t0, &v0, &tx0, &t1, &v1, &tx1, &t2, &v2, &tx2, &t3, &v3, &tx3);
            if (err == 9)
              {
              tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
              triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
              }
            else if (err == 12)
              {
              tria_uv.push_back(vec3<uint32_t>(v0 - 1, v1 - 1, v2 - 1));
              triangles.push_back(vec3<uint32_t>(t0 - 1, t1 - 1, t2 - 1));
              tria_uv.push_back(vec3<uint32_t>(v0 - 1, v2 - 1, v3 - 1));
              triangles.push_back(vec3<uint32_t>(t0 - 1, t2 - 1, t3 - 1));
              }
            else
              {
              fclose(f);
              return false;
              }

            }
          }
        }
      }
    }
  fclose(f);
  if (!tria_uv.empty() && (triangles.size() != tria_uv.size()))
    return false;
  if (!tria_uv.empty())
    {
    for (auto& t : tex)
      {
      t[1] = 1.f - t[1];
      }
    uv.reserve(triangles.size());
    for (auto t : tria_uv)
      {
      uv.emplace_back(tex[t[0]], tex[t[1]], tex[t[2]]);
      }
    }

  if (!mtl_filename.empty())
    {
    if (!file_exists(mtl_filename))
      {
      mtl_filename = get_folder(filename) + mtl_filename;
      }
    std::string texture_filename;
    if (read_texture_filename_from_mtl(texture_filename, mtl_filename.c_str()))
      {
      if (!file_exists(texture_filename))
        {
        texture_filename = get_folder(mtl_filename) + texture_filename;
        }
      int w, h, nr_of_channels;
      unsigned char* im = stbi_load(texture_filename.c_str(), &w, &h, &nr_of_channels, 4);
      if (im)
        {
        texture = jtk::span_to_image(w, h, w, (const uint32_t*)im);
        stbi_image_free(im);
        }
      }
    }

  return true;
  }

bool write_obj(const char* filename, const std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<float>>& normals, const std::vector<uint32_t>& clrs, const std::vector<jtk::vec3<uint32_t>>& triangles, const std::vector<jtk::vec3<jtk::vec2<float>>>& uv, const jtk::image<uint32_t>& texture)
  {
  using namespace jtk;
  std::string objpath(filename);

  std::string mtlfilename;

  if (texture.width() > 0 && texture.height() > 0)
    {

    std::string fn = get_filename(objpath);
    std::string folder = get_folder(objpath);
    std::string imagefilename = fn + ".png";
    mtlfilename = fn + ".mtl";

    std::string imagepath = folder.empty() ? imagefilename : folder + imagefilename;
    std::string mtlpath = folder.empty() ? mtlfilename : folder + mtlfilename;

    std::ofstream mtlfile;

    mtlfile.open(mtlpath.c_str(), std::ios::out);
    if (mtlfile.fail())
      return false;
    if (mtlfile.bad())
      return false;
    if (mtlfile.is_open())
      {
      mtlfile << "newmtl material_0" << std::endl;
      mtlfile << "Ka 0.0000 0.0000 0.0000" << std::endl;
      mtlfile << "Kd 0.8000 0.8000 0.8000" << std::endl;
      mtlfile << "Ks 1.0000 1.0000 1.0000" << std::endl;
      mtlfile << "Tf 0.0000 0.0000 0.0000" << std::endl;
      mtlfile << "d 1.0000" << std::endl;
      mtlfile << "Ns 0" << std::endl;
      mtlfile << "map_Kd " << imagefilename << std::endl;
      mtlfile.clear();
      }
    mtlfile.close();

    if (!stbi_write_png(imagepath.c_str(), texture.width(), texture.height(), 4, (void*)texture.data(), texture.stride() * 4))
      return false;
    }

  std::ofstream outputfile;
  outputfile.open(objpath.c_str(), std::ios::out);
  if (outputfile.fail())
    return false;
  if (outputfile.bad())
    return false;
  if (outputfile.is_open())
    {
    outputfile.clear();
    if (!mtlfilename.empty())
      {
      outputfile << "mtllib " << mtlfilename << std::endl;
      outputfile << "usemtl material_0" << std::endl;
      }
    for (size_t i = 0; i < vertices.size(); ++i)
      {
      outputfile << "v " << vertices[i][0] << " " << vertices[i][1] << " " << vertices[i][2] << std::endl;
      }
    for (size_t i = 0; i < normals.size(); ++i)
      {
      outputfile << "vn " << normals[i][0] << " " << normals[i][1] << " " << normals[i][2] << std::endl;
      }
    for (size_t i = 0; i < uv.size(); ++i)
      {
      outputfile << "vt " << uv[i][0][0] << " " << 1.f-uv[i][0][1] << std::endl;
      outputfile << "vt " << uv[i][1][0] << " " << 1.f-uv[i][1][1] << std::endl;
      outputfile << "vt " << uv[i][2][0] << " " << 1.f-uv[i][2][1] << std::endl;
      }
    for (size_t i = 0; i < triangles.size(); ++i)
      {
      if (uv.empty())
        {
        outputfile << "f " << triangles[i][0] + 1 << " " << triangles[i][1] + 1 << " " << triangles[i][2] + 1 << std::endl;
        }
      else
        {
        outputfile << "f " << triangles[i][0] + 1 << "/" << (i * 3) + 1 << " ";
        outputfile << triangles[i][1] + 1 << "/" << (i * 3 + 1) + 1 << " ";
        outputfile << triangles[i][2] + 1 << "/" << (i * 3 + 2) + 1 << std::endl;
        }
      }
    }
  outputfile.close();

  return true;
  }


bool read_pts(const char* filename, std::vector<jtk::vec3<float>>& vertices, std::vector<int>& intensity, std::vector<uint32_t>& clrs)
  {
  vertices.clear();
  intensity.clear();
  clrs.clear();

  FILE* f = nullptr;
  f = fopen(filename, "r");
  if (!f)
    return false;

  char buffer[1024];
  uint64_t nr_of_points = 0;

  while (fgets(buffer, 1024, f) != nullptr)
    {
    int first_non_whitespace_index = 0;
    while (first_non_whitespace_index < 250 && (buffer[first_non_whitespace_index] == ' ' || buffer[first_non_whitespace_index] == '\t'))
      ++first_non_whitespace_index;

    float x, y, z;
    int in;
    uint32_t r, g, b;    

    auto err = sscanf(buffer + first_non_whitespace_index, "%f %f %f %d %d %d %d\n", &x, &y, &z, &in, &r, &g, &b);
    if (err == 7)
      {
      vertices.emplace_back(x, y, z);
      intensity.push_back(in);
      uint32_t clr = 0xff000000 | (b << 16) | (g << 8) | r;
      clrs.push_back(clr);
      }
    else
      {
      err = sscanf(buffer + first_non_whitespace_index, "%f %f %f\n", &x, &y, &z);
      if (err == 3)
        {
        vertices.emplace_back(x, y, z);
        }
      else
        {
        err = sscanf(buffer + first_non_whitespace_index, "%d\n", &nr_of_points);
        if (err != 1)
          {
          std::cout << "Invalid pts file\n";
          return false;
          }
        }
      }
    }

  if (nr_of_points && vertices.size() != nr_of_points)
    {
    std::cout << "Invalid pts file\n";
    return false;
    }

  return true;
  }


bool write_pts(const char* filename, const std::vector<jtk::vec3<float>>& vertices, const std::vector<int>& intensity, const std::vector<uint32_t>& clrs)
  {
  std::ofstream out;

  out.open(filename, std::ios::out);
  if (out.fail())
    return false;
  if (out.bad())
    return false;
  if (out.is_open())
    {
    out << vertices.size() << std::endl;
    for (size_t i = 0; i < vertices.size(); ++i)
      {
      float x = vertices[i][0];
      float y = vertices[i][1];
      float z = vertices[i][2];
      int in = 0;
      if (intensity.size() == vertices.size())
        in = intensity[i];
      uint32_t clr = 0;
      if (clrs.size() == vertices.size())
        clr = clrs[i];
      uint32_t r = clr & 255;
      uint32_t g = (clr >> 8) & 255;
      uint32_t b = (clr >> 16) & 255;
      out << x << " " << y << " " << z << " " << in << " " << r << " " << g << " " << b << std::endl;
      }
    out.close();
    }
  return true;
  }

bool read_xyz(const char* filename, std::vector<jtk::vec3<float>>& vertices)
  {
  vertices.clear();
  std::ifstream in(filename, std::ios::in);

  if (in.fail())
    return false;
  if (in.bad())
    return false;
  if (in.is_open())
    {
    while (!in.eof())
      {
      std::string line;
      std::getline(in, line);
      if (!line.empty())
        {
        std::stringstream str;
        str << line;
        float x, y, z;
        str >> x >> y >> z;
        vertices.emplace_back(x, y, z);
        }
      }
    in.close();
    }
  return true;
  }

bool write_xyz(const char* filename, const std::vector<jtk::vec3<float>>& vertices)
  {
  std::ofstream out;

  out.open(filename, std::ios::out);
  if (out.fail())
    return false;
  if (out.bad())
    return false;
  if (out.is_open())
    {
    out << vertices.size() << std::endl;
    for (size_t i = 0; i < vertices.size(); ++i)
      {
      float x = vertices[i][0];
      float y = vertices[i][1];
      float z = vertices[i][2];     
      out << x << " " << y << " " << z << std::endl;
      }
    out.close();
    }
  return true;
  }