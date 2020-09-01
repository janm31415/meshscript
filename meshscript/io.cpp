#include "io.h"
#include <string.h>

extern "C"
  {
#include "ply.h"
  }

/*

namespace
  {

  struct vert
    {
    float x, y, z, nx, ny, nz;
    unsigned char r, g, b;
    };

  struct tria
    {
    unsigned char nverts;
    int* verts;
    };

  PlyProperty vert_props[] = { 
  {"x", PLY_FLOAT, PLY_FLOAT, offsetof(vert, x), 0, 0, 0, 0},
  {"y", PLY_FLOAT, PLY_FLOAT, offsetof(vert, y), 0, 0, 0, 0},
  {"z", PLY_FLOAT, PLY_FLOAT, offsetof(vert, z), 0, 0, 0, 0},
  {"nx", PLY_FLOAT, PLY_FLOAT, offsetof(vert, nx), 0, 0, 0, 0},
  {"ny", PLY_FLOAT, PLY_FLOAT, offsetof(vert, ny), 0, 0, 0, 0},
  {"nz", PLY_FLOAT, PLY_FLOAT, offsetof(vert, nz), 0, 0, 0, 0},
  {"diffuse_red", PLY_UCHAR, PLY_UCHAR, offsetof(vert, r), 0, 0, 0, 0},
  {"diffuse_green", PLY_UCHAR, PLY_UCHAR, offsetof(vert, g), 0, 0, 0, 0},
  {"diffuse_blue", PLY_UCHAR, PLY_UCHAR, offsetof(vert, b), 0, 0, 0, 0},
  {"red", PLY_UCHAR, PLY_UCHAR, offsetof(vert, r), 0, 0, 0, 0},
  {"green", PLY_UCHAR, PLY_UCHAR, offsetof(vert, g), 0, 0, 0, 0},
  {"blue", PLY_UCHAR, PLY_UCHAR, offsetof(vert, b), 0, 0, 0, 0},
    };

  PlyProperty face_props[] = {
      {"vertex_indices", PLY_INT, PLY_INT, offsetof(tria, verts), 1, PLY_UCHAR, PLY_UCHAR, offsetof(tria, nverts)}
    };

  }

bool read_ply(const char* filename, std::vector<jtk::vec3<float>>& pts, std::vector<jtk::vec3<float>>& normals, std::vector<uint32_t>& clrs)
  {
  pts.clear();
  normals.clear();
  clrs.clear();

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
      return 0;
      }
    if (equal_strings("vertex", elem_name))
      {
      for (int p = 0; p < nr_props; ++p)
        {
        for (int vp = 0; vp < 12; ++vp)
          {
          if (strcmp(plist[p]->name, vert_props[vp].name) == 0)
            {
            ply_get_property(ply, elem_name, &vert_props[vp]);
            break;
            }
          }

        }

      pts.resize(num_elems);
      normals.resize(num_elems);
      clrs.resize(num_elems);

      for (j = 0; j < num_elems; j++)
        {
        vert v;
        ply_get_element(ply, (void *)&v);
        pts[j][0] = v.x;
        pts[j][1] = v.y;
        pts[j][2] = v.z;
        normals[j][0] = v.nx;
        normals[j][1] = v.ny;
        normals[j][2] = v.nz;
        uint32_t r = (uint32_t)v.r;
        uint32_t g = (uint32_t)v.g;
        uint32_t b = (uint32_t)v.b;
        clrs[j] = 0xff000000 | (b << 16) | (g << 8) | r;
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
  return true;
  }


bool read_ply(const char* filename, std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles, std::vector<uint32_t>& clrs)
  {
  vertices.clear();
  clrs.clear();

  bool has_color = false;

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
      return 0;
      }
    if (equal_strings("vertex", elem_name))
      {
      for (int p = 0; p < nr_props; ++p)
        {
        for (int vp = 0; vp < 12; ++vp)
          {
          if (strcmp(plist[p]->name, vert_props[vp].name) == 0)
            {
            if (vp >= 6 && vp < 12)
              has_color = true;
            ply_get_property(ply, elem_name, &vert_props[vp]);
            break;
            }
          }

        }

      vertices.resize(num_elems);
      if (has_color)
        clrs.resize(num_elems);

      for (j = 0; j < num_elems; j++)
        {
        vert v;
        ply_get_element(ply, (void *)&v);
        vertices[j][0] = v.x;
        vertices[j][1] = v.y;
        vertices[j][2] = v.z;
        if (has_color)
          {
          uint32_t r = (uint32_t)v.r;
          uint32_t g = (uint32_t)v.g;
          uint32_t b = (uint32_t)v.b;
          clrs[j] = 0xff000000 | (b << 16) | (g << 8) | r;
          }
        }
      }
    if (equal_strings("face", elem_name))
      {
      for (int p = 0; p < nr_props; ++p)
        {
        if (strcmp(plist[p]->name, face_props[0].name) == 0)
          {
          ply_get_property(ply, elem_name, &face_props[0]);
          }
        }

      triangles.resize(num_elems);

      for (j = 0; j < num_elems; j++)
        {
        tria t;
        ply_get_element(ply, (void *)&t);

        if (t.nverts == 3)
          {
          triangles[j][0] = t.verts[0];
          triangles[j][1] = t.verts[1];
          triangles[j][2] = t.verts[2];
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
  return true;
  }

*/
//////////////////////////////////////
//////////////////////////////////////

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