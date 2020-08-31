#include "io.h"
#include <string.h>

extern "C"
  {
#include "ply.h"
  }

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

  PlyProperty vert_props[] = { /* list of property information for a vertex */
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