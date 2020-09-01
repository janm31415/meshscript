#include "io.h"
#include <string.h>

#include "rply.h"

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

  if (nred > 0)
    clrs.resize(nred);

  p_red = (uint8_t*)clrs.data();
  p_green = p_red + 1;
  p_blue = p_red + 2;
  p_alpha = p_red + 3;

  uint32_t* p_tria_index = NULL;

  long ntriangles = ply_set_read_cb(ply, "face", "vertex_indices", read_face, (void*)(&p_tria_index), 0);

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