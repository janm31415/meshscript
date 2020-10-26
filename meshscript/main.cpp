#ifndef _SDL_main_h
#define _SDL_main_h
#endif
#include <SDL.h>

#include <stdio.h>
#include <string.h>

#include <libskiwi/libskiwi.h>

#include <libskiwi/runtime.h>

#include <libskiwi/types.h>

#include "view.h"
#include "jet.h"
#include <thread>
#include <memory>
#include <mutex>
#include <condition_variable>
#include <sstream>
#include <iostream>

#include <jtk/concurrency.h>
#include <jtk/fitting.h>
#include <jtk/geometry.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

using namespace jtk;

constexpr double pi = 3.141592653589793238462643383279;

view* g_view = nullptr;

bool quit = false;

void scm_exit()
  {
  quit = true;
  }

void scm_hide_view()
  {
  g_view->hide();
  }

void scm_show_view()
  {
  g_view->show();
  }

int64_t load_mesh(const char* filename)
  {
  int64_t id = g_view->load_mesh_from_file(filename);
  return id;
  }

int64_t scm_mesh_to_pointcloud(int64_t id)
  {
  int64_t id_out = g_view->mesh_to_pointcloud((uint32_t)id);
  return id_out;
  }

int64_t load_morphable_model(const char* filename)
  {
  int64_t id = g_view->load_morphable_model_from_file(filename);
  return id;
  }

int64_t load_pc(const char* filename)
  {
  int64_t id = g_view->load_pc_from_file(filename);
  return id;
  }

/*
Input can be a vector of size 16 in column major format,
or a list of lists in row major format like ((1 0 0 13) (0 1 0 12) (0 0 1 15) (0 0 0 1))
which can be read from a csv file.
*/
void scm_set_coordinate_system(int64_t id, uint64_t scheme_variable_64)
  {
  skiwi::scm_type scheme_variable(scheme_variable_64);
  if (scheme_variable.is_vector())
    {
    auto v = scheme_variable.get_vector();
    if (v.size() != 16)
      {
      std::cout << "error: cs-set!: second input parameter should be a vector of size 16\n";
      std::cout << "                current size of vector is " << v.size() << "\n";
      return;
      }
    float4x4 m;
    try
      {
      for (int i = 0; i < 16; ++i)
        {
        m[i] = (float)v[i].get_number();
        }
      g_view->set_coordinate_system((uint32_t)id, m);
      }
    catch (std::runtime_error e)
      {
      std::cout << e.what() << "\n";
      }
    }
  else if (scheme_variable.is_pair())
    {
    try
      {
      std::vector<skiwi::scm_type> rows(4);
      rows[0] = scheme_variable.get_pair().first;
      rows[1] = scheme_variable.get_pair().second.get_pair().first;
      rows[2] = scheme_variable.get_pair().second.get_pair().second.get_pair().first;
      rows[3] = scheme_variable.get_pair().second.get_pair().second.get_pair().second.get_pair().first;
      float4x4 m;
      for (int i = 0; i < 4; ++i)
        {
        const auto& r = rows[i];
        auto v0 = r.get_pair().first;
        auto v1 = r.get_pair().second.get_pair().first;
        auto v2 = r.get_pair().second.get_pair().second.get_pair().first;
        auto v3 = r.get_pair().second.get_pair().second.get_pair().second.get_pair().first;
        m[i] = (float)v0.get_number();
        m[i + 4] = (float)v1.get_number();
        m[i + 8] = (float)v2.get_number();
        m[i + 12] = (float)v3.get_number();
        }
      g_view->set_coordinate_system((uint32_t)id, m);
      }
    catch (std::runtime_error e)
      {
      std::cout << e.what() << "\n";
      }
    }
  else
    {
    std::cout << "error: cs-set!: invalid input type\n";
    }
  }

/*
Input can be a vector of size 16 in column major format,
or a list of lists in row major format like ((1 0 0 13) (0 1 0 12) (0 0 1 15) (0 0 0 1))
which can be read from a csv file.
*/
void scm_set_view_coordinate_system(uint64_t scheme_variable_64)
  {
  skiwi::scm_type scheme_variable(scheme_variable_64);
  if (scheme_variable.is_vector())
    {
    auto v = scheme_variable.get_vector();
    if (v.size() != 16)
      {
      std::cout << "error: view-cs-set!: second input parameter should be a vector of size 16\n";
      std::cout << "                     current size of vector is " << v.size() << "\n";
      return;
      }
    float4x4 m;
    try
      {
      for (int i = 0; i < 16; ++i)
        {
        m[i] = (float)v[i].get_number();
        }
      g_view->set_coordinate_system(m);
      }
    catch (std::runtime_error e)
      {
      std::cout << e.what() << "\n";
      }
    }
  else if (scheme_variable.is_pair())
    {
    try
      {
      std::vector<skiwi::scm_type> rows(4);
      rows[0] = scheme_variable.get_pair().first;
      rows[1] = scheme_variable.get_pair().second.get_pair().first;
      rows[2] = scheme_variable.get_pair().second.get_pair().second.get_pair().first;
      rows[3] = scheme_variable.get_pair().second.get_pair().second.get_pair().second.get_pair().first;
      float4x4 m;
      for (int i = 0; i < 4; ++i)
        {
        const auto& r = rows[i];
        auto v0 = r.get_pair().first;
        auto v1 = r.get_pair().second.get_pair().first;
        auto v2 = r.get_pair().second.get_pair().second.get_pair().first;
        auto v3 = r.get_pair().second.get_pair().second.get_pair().second.get_pair().first;
        m[i] = (float)v0.get_number();
        m[i + 4] = (float)v1.get_number();
        m[i + 8] = (float)v2.get_number();
        m[i + 12] = (float)v3.get_number();
        }
      g_view->set_coordinate_system(m);
      }
    catch (std::runtime_error e)
      {
      std::cout << e.what() << "\n";
      }
    }
  else
    {
    std::cout << "error: view-cs-set!: invalid input type\n";
    }
  }

uint64_t scm_get_coordinate_system(int64_t id)
  {
  using namespace skiwi;
  auto cs = g_view->get_coordinate_system((uint32_t)id);
  std::vector<scm_type> lst;
  for (int r = 0; r < 4; ++r)
    {
    std::vector<scm_type> row;
    for (int c = 0; c < 4; ++c)
      {
      row.push_back(make_flonum(cs[r + 4 * c]));
      }
    lst.push_back(make_list(row));
    }
  return make_list(lst);
  }

uint64_t scm_icp(int64_t id1, int64_t id2, uint64_t inlier_distance)
  {
  using namespace skiwi;
  scm_type inlier_dist(inlier_distance);
  double inlierd = inlier_dist.get_number();
  auto cs = g_view->icp((uint32_t)id1, (uint32_t)id2, inlierd);
  std::vector<scm_type> lst;
  for (int r = 0; r < 4; ++r)
    {
    std::vector<scm_type> row;
    for (int c = 0; c < 4; ++c)
      {
      row.push_back(make_flonum(cs[r + 4 * c]));
      }
    lst.push_back(make_list(row));
    }
  return make_list(lst);
  }

uint64_t scm_distance_map(int64_t id1, int64_t id2, bool sign)
  {
  using namespace skiwi;
  auto lst = g_view->distance_map((uint32_t)id1, (uint32_t)id2, sign);
  std::vector<scm_type> outlist;
  for (const auto& v : lst)
    {
    outlist.push_back(make_flonum((double)v));
    }
  return make_list(outlist);
  }

uint64_t scm_get_view_coordinate_system()
  {
  using namespace skiwi;
  auto cs = g_view->get_coordinate_system();
  std::vector<scm_type> lst;
  for (int r = 0; r < 4; ++r)
    {
    std::vector<scm_type> row;
    for (int c = 0; c < 4; ++c)
      {
      row.push_back(make_flonum(cs[r + 4 * c]));
      }
    lst.push_back(make_list(row));
    }
  return make_list(lst);
  }

namespace
  {

  jtk::combinable<void*> thread_specific_scheme_context;

  skiwi::skiwi_compiled_function_ptr marching_cubes_fun = nullptr;

  double marching_cubes_distance_fun(double x, double y, double z)
    {
    bool exists;
    void*& local_context = thread_specific_scheme_context.local(exists);
    if (!exists)
      local_context = skiwi::skiwi_clone_context(skiwi::skiwi_get_context());
    skiwi::scm_type res = skiwi::skiwi_run_raw(marching_cubes_fun, local_context, x, y, z);
    return res.get_number();
    }

  }

int64_t scm_marching_cubes(uint64_t bb64, uint64_t dim64, uint64_t iso64, uint64_t fun64)
  {
  skiwi::scm_type bb(bb64);
  skiwi::scm_type dim(dim64);
  skiwi::scm_type iso(iso64);
  skiwi::scm_type fun(fun64);

  skiwi::save_compiler_data();

  jtk::boundingbox3d<float> bounding;
  int64_t width, height, depth;
  double isovalue;
  if (bb.is_pair() && dim.is_pair() && (iso.is_fixnum() || iso.is_flonum()) && fun.is_closure())
    {
    try
      {
      std::vector<skiwi::scm_type> bb_(3);
      bb_[0] = bb.get_pair().first;
      bb_[1] = bb.get_pair().second.get_pair().first;
      bb_[2] = bb.get_pair().second.get_pair().second.get_pair().first;

      for (int i = 0; i < 3; ++i)
        {
        const auto& minmax = bb_[i];
        auto mi = minmax.get_pair().first;
        auto ma = minmax.get_pair().second.get_pair().first;
        bounding.min[i] = (float)mi.get_number();
        bounding.max[i] = (float)ma.get_number();
        }

      width = dim.get_pair().first.get_fixnum();
      height = dim.get_pair().second.get_pair().first.get_fixnum();
      depth = dim.get_pair().second.get_pair().second.get_pair().first.get_fixnum();
      isovalue = iso.get_number();

      std::string script = std::string("(c-input \"(double mc_x, double mc_y, double mc_z)\") (") + fun.get_closure_name() + std::string(" mc_x mc_y mc_z)");

      marching_cubes_fun = skiwi::skiwi_compile(script);

      int64_t id = g_view->marching_cubes(bounding, width, height, depth, (float)isovalue, &marching_cubes_distance_fun);

      thread_specific_scheme_context.combine_each([](void* ctxt)
        {
        skiwi::skiwi_destroy_clone_context(ctxt);
        });
      thread_specific_scheme_context.clear();
      skiwi::restore_compiler_data();
      return id;
      }
    catch (std::runtime_error e)
      {
      std::cout << e.what() << "\n";
      }
    }
  else
    {
    std::cout << "error: marching-cubes: invalid input type\n";
    }
  skiwi::restore_compiler_data();
  return -1;
  }

int64_t make_mesh(uint64_t scm_vertices_64, uint64_t scm_triangles_64)
  {
  skiwi::scm_type scm_vertices(scm_vertices_64);
  skiwi::scm_type scm_triangles(scm_triangles_64);

  try
    {
    std::vector<vec3<float>> vertices;

    {
    auto vert = scm_vertices.get_list();
    vertices.reserve(vert.size());
    for (auto& v : vert)
      {
      auto vertex = v.get_list();
      if (vertex.size() != 3)
        throw std::runtime_error("error: make-mesh: invalid vertex size (should have 3 float values)");
      vertices.emplace_back((float)vertex[0].get_number(), (float)vertex[1].get_number(), (float)vertex[2].get_number());
      }
    }

    std::vector<vec3<uint32_t>> triangles;
    {
    auto tria = scm_triangles.get_list();
    triangles.reserve(tria.size());
    for (auto& t : tria)
      {
      auto triangle = t.get_list();
      if (triangle.size() != 3)
        throw std::runtime_error("error: make-mesh: invalid triangle size (should have 3 int values)");
      triangles.emplace_back((uint32_t)triangle[0].get_fixnum(), (uint32_t)triangle[1].get_fixnum(), (uint32_t)triangle[2].get_fixnum());
      }
    }
    int64_t id = g_view->load_mesh(vertices, triangles);
    return id;
    }
  catch (std::runtime_error e)
    {
    std::cout << e.what() << "\n";
    }
  return -1;
  }

bool vertices_to_csv(int64_t id, const char* filename)
  {
  bool res = g_view->vertices_to_csv(id, filename);
  return res;
  }

bool triangles_to_csv(int64_t id, const char* filename)
  {
  bool res = g_view->triangles_to_csv(id, filename);
  return res;
  }

void show(int64_t id)
  {
  g_view->show(id);
  }

void hide(int64_t id)
  {
  g_view->hide(id);
  }

void set_matcap(int64_t id, int64_t clr_id)
  {
  g_view->set_matcap(id, clr_id);
  }

void scm_unzoom()
  {
  g_view->unzoom();
  }

uint64_t scm_jet(uint64_t mask64)
  {
  using namespace skiwi;
  skiwi::scm_type mask(mask64);
  try
    {
    auto m = mask.get_list();
    std::vector<scm_type> colors;
    colors.reserve(m.size());
    for (auto value : m)
      {
      float gray = (float)value.get_number();
      if (gray < 0.f)
        gray = 0.f;
      if (gray > 1.f)
        gray = 1.f;

      //gray *= 2.f;
      //gray -= 1.f;

      float r, g, b;
      jet_cold_to_hot(r, g, b, gray);
      std::swap(g, b);
      //r = jet_red(gray);
      //b = jet_green(gray); // swap blue and green
      //g = jet_blue(gray);
      int64_t r2 = (int64_t)(r*255.f);
      int64_t g2 = (int64_t)(g*255.f);
      int64_t b2 = (int64_t)(b*255.f);
      std::vector<scm_type> clr;
      clr.push_back(make_fixnum(r2));
      clr.push_back(make_fixnum(g2));
      clr.push_back(make_fixnum(b2));
      colors.push_back(make_list(clr));
      }
    return make_list(colors);
    }
  catch (std::runtime_error e)
    {
    std::cout << "error: jet: " << e.what() << "\n";
    }
  return make_undefined();
  }

void scm_set_vertex_colors(int64_t id, uint64_t scm_colors_64)
  {
  skiwi::scm_type scm_colors(scm_colors_64);
  std::vector<vec3<uint8_t>> colors;
  try
    {
    auto clrlst = scm_colors.get_list();
    colors.reserve(clrlst.size());
    for (auto& clr : clrlst)
      {
      auto c = clr.get_list();
      if (c.size() != 3)
        throw std::runtime_error("error: set-vertex-colors: invalid color size (should have 3 int values)");
      colors.emplace_back((uint8_t)c[0].get_fixnum(), (uint8_t)c[1].get_fixnum(), (uint8_t)c[2].get_fixnum());
      }
    g_view->set_vertex_colors((uint32_t)id, colors);
    }
  catch (std::runtime_error e)
    {
    std::cout << "error: set-vertex-colors: " << e.what() << "\n";
    }
  }

void scm_rotate(int64_t id, uint64_t x_axis_64, uint64_t y_axis_64, uint64_t z_axis_64)
  {
  skiwi::scm_type x_axis(x_axis_64);
  skiwi::scm_type y_axis(y_axis_64);
  skiwi::scm_type z_axis(z_axis_64);
  double x = x_axis.get_number();
  double y = y_axis.get_number();
  double z = z_axis.get_number();
  x *= pi / 180.0;
  y *= pi / 180.0;
  z *= pi / 180.0;
  auto rot_x = jtk::make_rotation(vec3<float>(0, 0, 0), vec3<float>(1, 0, 0), (float)x);
  auto rot_y = jtk::make_rotation(vec3<float>(0, 0, 0), vec3<float>(0, 1, 0), (float)y);
  auto rot_z = jtk::make_rotation(vec3<float>(0, 0, 0), vec3<float>(0, 0, 1), (float)z);

  auto rot = jtk::matrix_matrix_multiply(jtk::matrix_matrix_multiply(rot_x, rot_y), rot_z);
  g_view->premultiply_coordinate_system((uint32_t)id, rot);
  }

void scm_cs_premultiply(int64_t id, uint64_t scheme_variable_64)
  {
  skiwi::scm_type scheme_variable(scheme_variable_64);
  if (scheme_variable.is_vector())
    {
    auto v = scheme_variable.get_vector();
    if (v.size() != 16)
      {
      std::cout << "error: cs-premultiply!: second input parameter should be a vector of size 16\n";
      std::cout << "                        current size of vector is " << v.size() << "\n";
      return;
      }
    float4x4 m;
    try
      {
      for (int i = 0; i < 16; ++i)
        {
        m[i] = (float)v[i].get_number();
        }
      g_view->premultiply_coordinate_system((uint32_t)id, m);
      }
    catch (std::runtime_error e)
      {
      std::cout << e.what() << "\n";
      }
    }
  else if (scheme_variable.is_pair())
    {
    try
      {
      std::vector<skiwi::scm_type> rows(4);
      rows[0] = scheme_variable.get_pair().first;
      rows[1] = scheme_variable.get_pair().second.get_pair().first;
      rows[2] = scheme_variable.get_pair().second.get_pair().second.get_pair().first;
      rows[3] = scheme_variable.get_pair().second.get_pair().second.get_pair().second.get_pair().first;
      float4x4 m;
      for (int i = 0; i < 4; ++i)
        {
        const auto& r = rows[i];
        auto v0 = r.get_pair().first;
        auto v1 = r.get_pair().second.get_pair().first;
        auto v2 = r.get_pair().second.get_pair().second.get_pair().first;
        auto v3 = r.get_pair().second.get_pair().second.get_pair().second.get_pair().first;
        m[i] = (float)v0.get_number();
        m[i + 4] = (float)v1.get_number();
        m[i + 8] = (float)v2.get_number();
        m[i + 12] = (float)v3.get_number();
        }
      g_view->premultiply_coordinate_system((uint32_t)id, m);
      }
    catch (std::runtime_error e)
      {
      std::cout << e.what() << "\n";
      }
    }
  else
    {
    std::cout << "error: cs-premultiply!: invalid input type\n";
    }
  }

void scm_translate(int64_t id, uint64_t x_axis_64, uint64_t y_axis_64, uint64_t z_axis_64)
  {
  skiwi::scm_type x_axis(x_axis_64);
  skiwi::scm_type y_axis(y_axis_64);
  skiwi::scm_type z_axis(z_axis_64);

  double x = x_axis.get_number();
  double y = y_axis.get_number();
  double z = z_axis.get_number();

  auto t = jtk::make_translation((float)x, (float)y, (float)z);
  g_view->premultiply_coordinate_system((uint32_t)id, t);
  }

void scm_set_shading(bool b)
  {
  g_view->set_shading(b);
  }

void scm_set_edges(bool b)
  {
  g_view->set_edges(b);
  }

void scm_set_textured(bool b)
  {
  g_view->set_textured(b);
  }

void scm_set_shadow(bool b)
  {
  g_view->set_shadow(b);
  }

void scm_set_wireframe(bool b)
  {
  g_view->set_wireframe(b);
  }

void scm_set_one_bit(bool b)
  {
  g_view->set_one_bit(b);
  }

void scm_set_image_size(uint64_t w64, uint64_t h64)
  {
  skiwi::scm_type w(w64);
  skiwi::scm_type h(h64);
  int width = (int)w.get_number();
  int height = (int)h.get_number();
  g_view->set_image_size(width, height);
  }

void scm_export_image(const char* filename)
  {
  auto im = g_view->get_image();
  stbi_write_png(filename, im.width(), im.height(), 4, (void*)im.data(), im.stride() * 4);
  }

void scm_set_bg_color(int64_t r, int64_t g, int64_t b)
  {
  g_view->set_bg_color((uint8_t)r, (uint8_t)g, (uint8_t)b);
  }

uint64_t scm_get_position(uint64_t x64, uint64_t y64)
  {
  using namespace skiwi;
  skiwi::scm_type x(x64);
  skiwi::scm_type y(y64);
  auto pos = g_view->get_world_position((int)x.get_number(), (int)y.get_number());
  std::vector<scm_type> coord;
  coord.push_back(make_flonum(pos[0]));
  coord.push_back(make_flonum(pos[1]));
  coord.push_back(make_flonum(pos[2]));
  return make_list(coord);
  }

uint64_t scm_get_index(uint64_t x64, uint64_t y64)
  {
  using namespace skiwi;
  skiwi::scm_type x(x64);
  skiwi::scm_type y(y64);
  auto idx = g_view->get_index((int)x.get_number(), (int)y.get_number());
  return make_fixnum(idx);
  }

uint64_t scm_get_id(uint64_t x64, uint64_t y64)
  {
  using namespace skiwi;
  skiwi::scm_type x(x64);
  skiwi::scm_type y(y64);
  auto id = g_view->get_id((int)x.get_number(), (int)y.get_number());
  return make_fixnum(id);
  }

void scm_force_redraw()
  {
  g_view->force_redraw();
  }

int64_t mm_coeff_size(int64_t id)
  {
  return g_view->mm_coeff_size((uint32_t)id);
  }

int64_t mm_shape_size(int64_t id)
  {
  return g_view->mm_shape_size((uint32_t)id);
  }

double mm_sigma(int64_t id, int64_t idx)
  {
  return g_view->mm_sigma((uint32_t)id, idx);
  }

uint64_t mm_coeff(int64_t id)
  {
  using namespace skiwi;
  std::vector<float> coeff = g_view->mm_coeff((uint32_t)id);
  std::vector<scm_type> coefflist;
  for (const auto& v : coeff)
    {
    coefflist.push_back(make_flonum((double)v));
    }
  return make_list(coefflist);
  }

uint64_t mm_basic_shape_coeff(int64_t id, int64_t shape_id)
  {
  using namespace skiwi;
  std::vector<float> coeff = g_view->mm_basic_shape_coeff((uint32_t)id, shape_id);
  std::vector<scm_type> coefflist;
  for (const auto& v : coeff)
    {
    coefflist.push_back(make_flonum((double)v));
    }
  return make_list(coefflist);
  }

void mm_coeff_set(int64_t id, uint64_t scm_coeff_64)
  {
  skiwi::scm_type scm_coeff(scm_coeff_64);
  std::vector<float> coeff;
  try
    {
    auto clrlst = scm_coeff.get_list();
    coeff.reserve(clrlst.size());
    for (auto& clr : clrlst)
      {
      double c = clr.get_number();
      coeff.push_back((float)c);
      }
    g_view->mm_coeff_set((uint32_t)id, coeff);
    }
  catch (std::runtime_error e)
    {
    std::cout << "error: morphable-model-coeff-set!: " << e.what() << "\n";
    }
  }

int64_t mm_color_coeff_size(int64_t id)
  {
  return g_view->mm_color_coeff_size((uint32_t)id);
  }

int64_t mm_color_shape_size(int64_t id)
  {
  return g_view->mm_color_shape_size((uint32_t)id);
  }

double mm_color_sigma(int64_t id, int64_t idx)
  {
  return g_view->mm_color_sigma((uint32_t)id, idx);
  }

uint64_t mm_color_coeff(int64_t id)
  {
  using namespace skiwi;
  std::vector<float> coeff = g_view->mm_color_coeff((uint32_t)id);
  std::vector<scm_type> coefflist;
  for (const auto& v : coeff)
    {
    coefflist.push_back(make_flonum((double)v));
    }
  return make_list(coefflist);
  }

uint64_t mm_color_basic_shape_coeff(int64_t id, int64_t shape_id)
  {
  using namespace skiwi;
  std::vector<float> coeff = g_view->mm_color_basic_shape_coeff((uint32_t)id, shape_id);
  std::vector<scm_type> coefflist;
  for (const auto& v : coeff)
    {
    coefflist.push_back(make_flonum((double)v));
    }
  return make_list(coefflist);
  }

void mm_color_coeff_set(int64_t id, uint64_t scm_coeff_64)
  {
  skiwi::scm_type scm_coeff(scm_coeff_64);
  std::vector<float> coeff;
  try
    {
    auto clrlst = scm_coeff.get_list();
    coeff.reserve(clrlst.size());
    for (auto& clr : clrlst)
      {
      double c = clr.get_number();
      coeff.push_back((float)c);
      }
    g_view->mm_color_coeff_set((uint32_t)id, coeff);
    }
  catch (std::runtime_error e)
    {
    std::cout << "error: morphable-model-color-coeff-set!: " << e.what() << "\n";
    }
  }

int64_t mm_to_mesh(int64_t id)
  {
  return g_view->mm_to_mesh((uint32_t)id);
  }

uint64_t scm_triangles(int64_t id)
  {
  using namespace skiwi;
  auto trias = g_view->triangles((uint32_t)id);
  std::vector<scm_type> trialist;
  for (const auto& tria : trias)
    {
    std::vector<scm_type> t;
    t.push_back(make_fixnum(tria[0]));
    t.push_back(make_fixnum(tria[1]));
    t.push_back(make_fixnum(tria[2]));
    trialist.push_back(make_list(t));
    }
  return make_list(trialist);
  }

uint64_t scm_vertices(int64_t id)
  {
  using namespace skiwi;
  auto verts = g_view->vertices((uint32_t)id);
  std::vector<scm_type> vertlist;
  for (const auto& vert : verts)
    {
    std::vector<scm_type> v;
    v.push_back(make_flonum(vert[0]));
    v.push_back(make_flonum(vert[1]));
    v.push_back(make_flonum(vert[2]));
    vertlist.push_back(make_list(v));
    }
  return make_list(vertlist);
  }

bool scm_write(int64_t id, const char* filename)
  {
  return g_view->write((uint32_t)id, filename);
  }

uint64_t scm_mesh_texture_to_vertexcolors(uint64_t id)
  {
  using namespace skiwi;
  std::vector<vec3<uint8_t>> colors = g_view->mesh_texture_to_vertexcolors((uint32_t)id);
  std::vector<scm_type> vertclrlist;
  for (const auto& clr : colors)
    {
    std::vector<scm_type> v;
    v.push_back(make_fixnum(clr[0]));
    v.push_back(make_fixnum(clr[1]));
    v.push_back(make_fixnum(clr[2]));
    vertclrlist.push_back(make_list(v));
    }
  return make_list(vertclrlist);
  }

int64_t scm_load_shape_predictor(const char* filename)
  {
  return g_view->load_shape_predictor(filename);
  }

void scm_set_show_ear_right_detector(bool b)
  {
  g_view->set_show_ear_right_detector(b);
  }

void scm_set_show_ear_left_detector(bool b)
  {
  g_view->set_show_ear_left_detector(b);
  }

void scm_set_show_face_detector(bool b)
  {
  g_view->set_show_face_detector(b);
  }

uint64_t scm_right_ear_detect()
  {
  using namespace skiwi;
  std::vector<rect> rectangles = g_view->ear_detect(true);
  std::vector<skiwi::scm_type> vec;
  vec.reserve(rectangles.size());
  for (const auto& pr : rectangles)
    {
    std::vector<skiwi::scm_type> r;
    r.push_back(make_fixnum(pr.x));
    r.push_back(make_fixnum(pr.y));
    r.push_back(make_fixnum(pr.w));
    r.push_back(make_fixnum(pr.h));
    vec.push_back(make_list(r));
    }
  return make_list(vec);
  }

uint64_t scm_left_ear_detect()
  {
  using namespace skiwi;
  std::vector<rect> rectangles = g_view->ear_detect(false);
  std::vector<skiwi::scm_type> vec;
  vec.reserve(rectangles.size());
  for (const auto& pr : rectangles)
    {
    std::vector<skiwi::scm_type> r;
    r.push_back(make_fixnum(pr.x));
    r.push_back(make_fixnum(pr.y));
    r.push_back(make_fixnum(pr.w));
    r.push_back(make_fixnum(pr.h));
    vec.push_back(make_list(r));
    }
  return make_list(vec);
  }

uint64_t scm_face_detect()
  {
  using namespace skiwi;
  std::vector<rect> rectangles = g_view->face_detect();
  std::vector<skiwi::scm_type> vec;
  vec.reserve(rectangles.size());
  for (const auto& pr : rectangles)
    {
    std::vector<skiwi::scm_type> r;
    r.push_back(make_fixnum(pr.x));
    r.push_back(make_fixnum(pr.y));
    r.push_back(make_fixnum(pr.w));
    r.push_back(make_fixnum(pr.h));
    vec.push_back(make_list(r));
    }
  return make_list(vec);
  }

void scm_shape_predictor_horizontal_flip_set(int64_t id, bool f)
  {
  g_view->shape_predictor_set_flip_horizontal((uint32_t)id, f);
  }

void scm_sp_link_to_face(int64_t id)
  {
  g_view->shape_predictor_link((uint32_t)id, sp::odl_facial);
  }

void scm_sp_link_to_ear_right(int64_t id)
  {
  g_view->shape_predictor_link((uint32_t)id, sp::odl_ear_right);
  }

void scm_sp_link_to_ear_left(int64_t id)
  {
  g_view->shape_predictor_link((uint32_t)id, sp::odl_ear_left);
  }

void scm_sp_link_remove(int64_t id)
  {
  g_view->shape_predictor_link((uint32_t)id, sp::odl_none);
  }

uint64_t scm_shape_predict(int64_t id, uint64_t rect64)
  {
  using namespace skiwi;
  skiwi::scm_type rectangles(rect64);
  if (rectangles.is_nil())
    return make_nil();
  if (!rectangles.is_pair())
    {
    std::cout << "error: shape-predict: I expect a list of the form (x y w h) as input.\n";
    return make_nil();
    }
  try
    {
    std::vector<skiwi::scm_type> list_of_rectangles;
    if (rectangles.get_pair().first.is_pair())
      list_of_rectangles = rectangles.get_list();
    else
      list_of_rectangles.push_back(rectangles);

    std::vector<std::pair<long, long>> landmarks;
    for (const auto& rectangle : list_of_rectangles)
      {
      std::vector<skiwi::scm_type> rect_values(4);
      rect_values[0] = rectangle.get_pair().first;
      rect_values[1] = rectangle.get_pair().second.get_pair().first;
      rect_values[2] = rectangle.get_pair().second.get_pair().second.get_pair().first;
      rect_values[3] = rectangle.get_pair().second.get_pair().second.get_pair().second.get_pair().first;
      std::vector<int> rect_values_int(4);
      for (int i = 0; i < 4; ++i)
        {
        rect_values_int[i] = rect_values[i].get_fixnum();
        }
      rect r;
      r.x = rect_values_int[0];
      r.y = rect_values_int[1];
      r.w = rect_values_int[2];
      r.h = rect_values_int[3];
      std::vector<std::pair<long, long>> lms = g_view->shape_predict((uint32_t)id, r);
      landmarks.insert(landmarks.end(), lms.begin(), lms.end());
      }
    std::vector<skiwi::scm_type> vec;
    vec.reserve(landmarks.size());
    for (const auto& pr : landmarks)
      vec.push_back(make_pair(make_fixnum(pr.first), make_pair(make_fixnum(pr.second), make_nil())));
    return make_list(vec);
    }
  catch (std::runtime_error e)
    {
    std::cout << "error: shape-predict: I expect a list of the form (x y w h) as input.\n";
    }

  return make_nil();
  }

uint64_t npoint_scm(uint64_t src64, uint64_t tgt64)
  {
  using namespace skiwi;
  skiwi::scm_type src(src64);
  skiwi::scm_type tgt(tgt64);
  try
    {
    auto src_list = src.get_list();
    auto tgt_list = tgt.get_list();
    jtk::mat source(src_list.size(), 3);
    jtk::mat target(tgt_list.size(), 3);
    if (source.rows() != target.rows())
      {
      std::cout << "error: npoint: source and target list have different length\n";
      return make_undefined();
      }
    uint64_t nr_rows = 0;
    for (int i = 0; i < source.rows(); ++i)
      {
      auto src_row = src_list[i].get_list();
      auto tgt_row = tgt_list[i].get_list();
      if (src_row.size() != 3 || tgt_row.size() != 3)
        {
        std::cout << "error: npoint: source and target list should contain lists of size 3 as elements\n";
        return make_undefined();
        }
      source(nr_rows, 0) = src_row[0].get_number();
      source(nr_rows, 1) = src_row[1].get_number();
      source(nr_rows, 2) = src_row[2].get_number();
      target(nr_rows, 0) = tgt_row[0].get_number();
      target(nr_rows, 1) = tgt_row[1].get_number();
      target(nr_rows, 2) = tgt_row[2].get_number();
      if (!std::isnan(source(nr_rows, 0))
        && !std::isnan(source(nr_rows, 1))
        && !std::isnan(source(nr_rows, 2))
        && !std::isnan(target(nr_rows, 0))
        && !std::isnan(target(nr_rows, 1))
        && !std::isnan(target(nr_rows, 2)))
      ++nr_rows;
      }
    source.resize(nr_rows, 3);
    target.resize(nr_rows, 3);
    auto m = jtk::npoint(source, target, false, true);
    std::vector<scm_type> m_row(4);
    std::vector<scm_type> rows;
    for (int j = 0; j < 4; ++j)
      {
      m_row[0] = make_flonum(m(j, 0));
      m_row[1] = make_flonum(m(j, 1));
      m_row[2] = make_flonum(m(j, 2));
      m_row[3] = make_flonum(m(j, 3));
      rows.push_back(make_list(m_row));
      }
    return make_list(rows);
    }
  catch (std::runtime_error e)
    {
    std::cout << "error: npoint: " << e.what() << "\n";
    }
  return make_undefined();
  }

void mm_fit_indices(int64_t mm_id, uint64_t indices64, uint64_t positions64)
  {
  skiwi::scm_type indices(indices64);
  skiwi::scm_type positions(positions64);
  try
    {
    auto indices_list = indices.get_list();
    std::vector<uint32_t> ind;
    for (auto& nr : indices_list)
      {
      ind.push_back((uint32_t)nr.get_number());
      }
    std::vector<vec3<float>> pos;
    auto pos_list = positions.get_list();
    for (auto& p : pos_list)
      {
      auto p_list = p.get_list();
      if (p_list.size() != 3)
        throw std::runtime_error("error: morphable-model-fit-indices!: invalid vertex size (should have 3 float values)");
      pos.emplace_back((float)p_list[0].get_number(), (float)p_list[1].get_number(), (float)p_list[2].get_number());
      }
    g_view->fit_mm_to_partial_positions((uint32_t)mm_id, ind, pos);
    }
  catch (std::runtime_error e)
    {
    std::cout << "error: morphable-model-fit-indices!: " << e.what() << "\n";
    }
  }

void mm_fit(int64_t mm_id, int64_t mesh_id)
  {  
  g_view->fit_mm((uint32_t)mm_id, (uint32_t)mesh_id);
  }

int64_t scm_poisson(int64_t id, int64_t depth)
  {
  return g_view->poisson((uint32_t)id, (uint32_t)depth);
  }

void scm_info(int64_t id)
  {
  return g_view->info((uint32_t)id);
  }

void scm_cs_apply(int64_t id)
  {
  return g_view->cs_apply((uint32_t)id);
  }

void* register_functions(void*)
  {
  using namespace skiwi;
  register_external_primitive("cs", (void*)&scm_get_coordinate_system, skiwi_scm, skiwi_int64, "(cs id) returns the coordinate system for the object with tag `id`.");
  register_external_primitive("cs-apply!", (void*)&scm_cs_apply, skiwi_void, skiwi_int64, "(cs-apply! id) transforms the vertices of object with tag `id` by its coordinate system, and sets its coordinate system to the world.");
  register_external_primitive("cs-rotate!", (void*)&scm_rotate, skiwi_void, skiwi_int64, skiwi_scm, skiwi_scm, skiwi_scm, "(cs-rotate! id x y z) rotates the object with tag `id` by `x` degrees over the x-axis, by `y` degrees over the y-axis, and by `z` degrees over the z-axis.");
  register_external_primitive("cs-set!", (void*)&scm_set_coordinate_system, skiwi_void, skiwi_int64, skiwi_scm, "(cs-set! id cs) sets a new coordinate system for the object with tag `id`. The coordinate system `cs` can be given as a vector of size 16 in column major format or as a list of lists in row major format.");
  register_external_primitive("cs-translate!", (void*)&scm_translate, skiwi_void, skiwi_int64, skiwi_scm, skiwi_scm, skiwi_scm, "(cs-translate! id x y z) translates the object with tag `id` by vector (x y z).");
  register_external_primitive("cs-premultiply!", (void*)&scm_cs_premultiply, skiwi_void, skiwi_int64, skiwi_scm, "(cs-premultiply! id cs) premultiplies the coordinate system of the object with tag `id` by the input coordinate system. The coordinate system `cs` can be given as a vector of size 16 in column major format or as a list of lists in row major format.");


  register_external_primitive("distance-map", (void*)&scm_distance_map, skiwi_scm, skiwi_int64, skiwi_int64, skiwi_bool, "(distance-map id1 id2 bool-signed) returns a list with values that represent the distance between objects with tag `id1` and `id2`. For each vertex of object `id1` there is exactly one distance in the list. The distance can be signed or unsigned, depending on the boolean value that is given to `bool-signed`.");

  register_external_primitive("ear-right-detect", (void*)&scm_right_ear_detect, skiwi_scm, "(ear-right-detect) runs the ear detector on the current view and returns a list of lists of the form ((x y w h) ...) where (x y w h) represents a rectangle containing the right ear starting in corner (x,y) and with sizes (w,h).");

  register_external_primitive("ear-left-detect", (void*)&scm_left_ear_detect, skiwi_scm, "(ear-left-detect) runs the ear detector on the current view and returns a list of lists of the form ((x y w h) ...) where (x y w h) represents a rectangle containing the left ear starting in corner (x,y) and with sizes (w,h).");

  register_external_primitive("face-detect", (void*)&scm_face_detect, skiwi_scm, "(face-detect) runs the face detector on the current view and returns a list of lists of the form ((x y w h) ...) where (x y w h) represents a rectangle containing the face starting in corner (x,y) and with sizes (w,h).");

  register_external_primitive("force-redraw", (void*)&scm_force_redraw, skiwi_void, "(force-redraw) redraws the canvas. This is useful if you want to use view-position in your script, as view-position uses the data of the last render of the view.");

  register_external_primitive("hide!", (void*)&hide, skiwi_void, skiwi_int64, "(hide! id) makes the object with tag `id` invisible.");
  register_external_primitive("icp", (void*)&scm_icp, skiwi_scm, skiwi_int64, skiwi_int64, skiwi_scm, "(icp id1 id2 inlier-distance) returns the result of the iterative closest point algorithm between objects with tag `id1` and `id2`. This result is always a 4x4 transformation matrix. The iterative closest point algorithm will only use correspondences between `id1` and `id2` if their distance is smaller than `inlier-distance`.");
  register_external_primitive("info", (void*)&scm_info, skiwi_void, skiwi_int64, "(info id) prints info on the object with tag `id`.");
  register_external_primitive("jet", (void*)&scm_jet, skiwi_scm, skiwi_scm, "(jet lst) takes a list `lst` of values between 0 and 1 and returns a list of lists with (r g b) values.");

  register_external_primitive("load-mesh", (void*)&load_mesh, skiwi_int64, skiwi_char_pointer, "(load-mesh \"stlfile.stl\") loads the STL file and returns an id. Similarly (load-mesh \"objfile.obj\") loads an OBJ file and returns the id. Other input mesh formats that are implemented are PLY and OFF.");
  register_external_primitive("load-morphable-model", (void*)&load_morphable_model, skiwi_int64, skiwi_char_pointer, "(load-morphable-model \"model2019_fullHead.h5\") loads morphable models following the hdf5 file format as used by the Basel Face Model project (https://faces.dmi.unibas.ch/bfm/bfm2019.html). The other file format that can be read is meshscripts own binary morphable model file format with extension SSM.");
  register_external_primitive("load-pointcloud", (void*)&load_pc, skiwi_int64, skiwi_char_pointer, "(load-pointcloud \"pointcloud.ply\") loads the PLY file as point cloud and returns an id. Other file formats are not yet supported.");
  register_external_primitive("load-shape-predictor", (void*)&scm_load_shape_predictor, skiwi_int64, skiwi_char_pointer, "(load-shape-predictor \"filename\") initializes the shape predictor with the data given by \"filename\" and returns the id. This is the dlib shape predictor (http://dlib.net). The 68 points facial landmarks predictor data can be downloaded from https://github.com/davisking/dlib-models");

  register_external_primitive("make-mesh", (void*)&make_mesh, skiwi_int64, skiwi_scm, skiwi_scm, "(make-mesh vertices triangles) creates the mesh with given `vertices` and `triangles`, and returns the id of the created object. `vertices` should be a list of lists of the form ((x y z) (x y z) ...) with x,y,z floating point values, and `triangles` should be a list of lists of the form ((a b c) (d e f) ...) with a,b... fixnums referring to the vertex indices.");
  register_external_primitive("marching-cubes", (void*)&scm_marching_cubes, skiwi_int64, skiwi_scm, skiwi_scm, skiwi_scm, skiwi_scm, "(marching-cubes bb dim isovalue fun) with `bb` representing the bounding box of the form ((min_x max_x) (min_y max_y) (min_z max_z)), `dim` representing the dimensions of the form (width height depth), `isovalue` a flonum representing the signed distance requested, and `fun` representing the distance functions as a lambda function accepting (x y z) values and returning a distance.");
  register_external_primitive("matcap-set!", (void*)&set_matcap, skiwi_void, skiwi_int64, skiwi_int64, "(matcap-set! id matcap-id) changes the matcap of the object with tag `id`. The matcap is given by its id matcap-id. Currently the matcaps in meshscript are hardcoded. There are 4 available matcaps with ids 0, 1, 2, 3.");

  register_external_primitive("mesh->pointcloud", (void*)&scm_mesh_to_pointcloud, skiwi_int64, skiwi_int64, "(mesh->pointcloud id) converts the mesh with tag `id` to a pointcloud.");
  register_external_primitive("mesh-texture->vertexcolors", (void*)&scm_mesh_texture_to_vertexcolors, skiwi_scm, skiwi_int64, "(mesh-texture->vertexcolors id) will return a list of lists of the form ((r g b) (r g b) ... ). Each vertex of the object with tag `id` has a corresponding (r g b) value. This (r g b) value is obtained from the texture of `id`, if available.");


  register_external_primitive("morphable-model-coefficients-size", (void*)&mm_coeff_size, skiwi_int64, skiwi_int64, "(morphable-model-coefficients-size mm_id) returns the number of coefficients for the morphable model with tag `mm_id`.");
  register_external_primitive("morphable-model-shape-size", (void*)&mm_shape_size, skiwi_int64, skiwi_int64, "(morphable-model-shape_size mm_id) returns the shape size for the morphable model with tag `mm_id`. This is equal to the number of rows in the U matrix, where a shape S is represented as S = mu + U*c, with mu the average shape, and c the coefficients vector.");
  register_external_primitive("morphable-model-sigma", (void*)&mm_sigma, skiwi_double, skiwi_int64, skiwi_int64, "(morphable-model-sigma mm_id idx) returns sigma for the morphable model with tag `mm_id` at coefficient index `idx`.");
  register_external_primitive("morphable-model-coefficients", (void*)&mm_coeff, skiwi_scm, skiwi_int64, "(morphable-model-coefficients mm_id) returns the list of coefficients for the morphable model with tag `mm_id`.");
  register_external_primitive("morphable-model-basic-shape-coefficients", (void*)&mm_basic_shape_coeff, skiwi_scm, skiwi_int64, skiwi_int64, "(morphable-model-basic-shape-coefficients mm_id idx) returns the list of coefficients of the `idx`-th shape that was used to generate this morphable model. Not all morphable models have this data. For instance the Basel shape model does not contain this data.");
  register_external_primitive("morphable-model-coefficients-set!", (void*)&mm_coeff_set, skiwi_void, skiwi_int64, skiwi_scm, "(morphable-model-coefficients-set! mm_id coeff) sets the list of coefficients for the morphable model with tag `mm_id`. Here `coeff` is a list of coefficient values, and its size should equal (morphable-model-coefficients-size mm_id).");
  register_external_primitive("morphable-model->mesh", (void*)&mm_to_mesh, skiwi_int64, skiwi_int64, "(morphable-model->mesh mm_id) converts the morphable model with tag `mm_id` to a mesh and returns the new id.");
  register_external_primitive("morphable-model-color-coefficients-size", (void*)&mm_color_coeff_size, skiwi_int64, skiwi_int64, "(morphable-model-color-coefficients-size mm_id) returns the number of color coefficients for the morphable model with tag `mm_id`.");
  register_external_primitive("morphable-model-color-shape-size", (void*)&mm_color_shape_size, skiwi_int64, skiwi_int64, "(morphable-model-color-shape_size mm_id) returns the shape size for the color part of the morphable model with tag `mm_id`. This is equal to the number of rows in the U_color matrix, where a the shape colors S_color are represented as S_color = mu_color + U_color*c, with mu_color the average shape colors, and c the color coefficients vector.");
  register_external_primitive("morphable-model-color-sigma", (void*)&mm_color_sigma, skiwi_double, skiwi_int64, skiwi_int64, "(morphable-model-color-sigma mm_id idx) returns sigma for the colors of the morphable model with tag `mm_id` at color coefficient index `idx`.");
  register_external_primitive("morphable-model-color-coefficients", (void*)&mm_color_coeff, skiwi_scm, skiwi_int64, "(morphable-model-color-coefficients mm_id) returns the list of color coefficients for the morphable model with tag `mm_id`.");
  register_external_primitive("morphable-model-color-basic-shape-coefficients", (void*)&mm_color_basic_shape_coeff, skiwi_scm, skiwi_int64, skiwi_int64, "(morphable-model-color-basic-shape-coefficients mm_id idx) returns the list of color coefficients of the `idx`-th color shape that was used to generate this morphable model. Not all morphable models have this data. For instance the Basel shape model does not contain this data.");
  register_external_primitive("morphable-model-color-coefficients-set!", (void*)&mm_color_coeff_set, skiwi_void, skiwi_int64, skiwi_scm, "(morphable-model-color-coefficients-set! mm_id coeff) sets the list of color coefficients for the morphable model with tag `mm_id`. Here `coeff` is a list of color coefficient values, and its size should equal (morphable-model-color-coefficients-size mm_id).");

  
  register_external_primitive("morphable-model-fit-indices!", (void*)&mm_fit_indices, skiwi_void, skiwi_int64, skiwi_scm, skiwi_scm, "(morphable-model-fit-indices! mm_id indices positions)");
  register_external_primitive("morphable-model-fit!", (void*)&mm_fit, skiwi_void, skiwi_int64, skiwi_int64, "(morphable-model-fit! mm_id mesh_id)");


  register_external_primitive("npoint", (void*)&npoint_scm, skiwi::skiwi_scm, skiwi::skiwi_scm, skiwi::skiwi_scm, "(npoint from to) computes the npoint-registration of the set of 3d points in `from` to the set of 3d points in `to`. The result is a 4x4 transformation matrix. Here `from` and `to` are lists of lists of the form ((x y z) (x y z) ...) and `from` and `to` should have the same amount of 3d points.");

  register_external_primitive("poisson", (void*)&scm_poisson, skiwi::skiwi_int64, skiwi::skiwi_int64, skiwi::skiwi_int64, "(poisson pc_id depth) applies Poisson surface reconstruction to the pointcloud with tag `pc_id`. This is the screened Poisson surface reconstruction algorithm by M. Kazhdan and H. Hoppe. You have to provide the parameter `depth` which represents the depth of the octree during Poisson surface reconstruction.");

  register_external_primitive("save", (void*)&scm_write, skiwi_bool, skiwi_int64, skiwi_char_pointer, "(save id \"file.ext\") writes the object with tag `id` to file. The filetype is determined by the extension that is given. You can export meshes to STL or PLY, pointclouds to PLY, morphable models to SSM."); // don't use write: gives naming conflict with slib

  register_external_primitive("shape-predict", (void*)&scm_shape_predict, skiwi_scm, skiwi_int64, skiwi_scm, "(shape-predict sp_id (x y w h)) or (shape-predict sp_id ((x y w h) ...)) runs the shape predictor with tag `sp_id` on the region defined by (x y w h) or on the regions defined by ((x y w h) ...) in the current view and returns the coordinates of the landmarks as a list of lists. The predictor should be initialized with load-shape-predictor.");
  register_external_primitive("shape-predictor-horizontal-flip-set!", (void*)&scm_shape_predictor_horizontal_flip_set, skiwi_void, skiwi_int64, skiwi_bool, "(shape-predictor-horizontal-flip-set! id #t/#f) toggles horizontal flipping of the shape predictor given by tag `id`.");
  register_external_primitive("shape-predictor-link-to-ear-left-detector", (void*)&scm_sp_link_to_ear_left, skiwi_void, skiwi_int64, "(shape-predictor-link-to-ear-left-detector id) links the shape predictor given by tag `id` to the ear left detector. The result of this operation is that the shape predictor is rendered automatically when the left ear detector's automatic rendering is on. You can turn on/off automatic rendering of the left ear detector with the command (view-ear-left-detector-set! #t/#f).");
  register_external_primitive("shape-predictor-link-to-ear-right-detector", (void*)&scm_sp_link_to_ear_right, skiwi_void, skiwi_int64, "(shape-predictor-link-to-ear-right-detector id) links the shape predictor given by tag `id` to the ear right detector. The result of this operation is that the shape predictor is rendered automatically when the right ear detector's automatic rendering is on. You can turn on/off automatic rendering of the right ear detector with the command (view-ear-right-detector-set! #t/#f).");
  register_external_primitive("shape-predictor-link-to-face-detector", (void*)&scm_sp_link_to_face, skiwi_void, skiwi_int64, "(shape-predictor-link-to-face-detector id) links the shape predictor given by tag `id` to the face detector. The result of this operation is that the shape predictor is rendered automatically when the face detector's automatic rendering is on. You can turn on/off automatic rendering of the face detector with the command (view-face-detector-set! #t/#f).");
  register_external_primitive("shape-predictor-unlink", (void*)&scm_sp_link_remove, skiwi_void, skiwi_int64, "(shape-predictor-unlink id) unlinks the shape predictor given by tag `id`, see shape-predictor-link-to-face-detector, shape-predictor-link-to-ear-right-detector, or shape-predictor-link-to-ear-left-detector.");

  register_external_primitive("show!", (void*)&show, skiwi_void, skiwi_int64, "(show! id) makes the object with tag `id` visible.");
  register_external_primitive("triangles", (void*)&scm_triangles, skiwi_scm, skiwi_int64, "(triangles id) returns the triangles of object with tag `id` as a list of lists of the form ((v0 v1 v2) (v3 v4 v4) ...) where each sublist (v0 v1 v2) contain the indices of the vertices that form a triangle. The actual vertex positions can be obtained with the command (vertices id).");
  register_external_primitive("triangles->csv", (void*)&triangles_to_csv, skiwi_bool, skiwi_int64, skiwi_char_pointer, "(triangles->csv id \"file.csv\") exports the triangles of the object with tag `id` to a csv file.");
  register_external_primitive("vertexcolors-set!", (void*)&scm_set_vertex_colors, skiwi_void, skiwi_int64, skiwi_scm, "(vertexcolors-set! id clrlst) sets vertex colors for the object with tag `id`. The vertex colors are given as a list of lists with (r g b) values.");
  register_external_primitive("vertices", (void*)&scm_vertices, skiwi_scm, skiwi_int64, "(vertices id) returns the vertices of object with tag `id` as a list of lists of the form ((x y z) (x y z) ...) where each sublist (x y z) is a 3d point representing the position of that vertex.");
  register_external_primitive("vertices->csv", (void*)&vertices_to_csv, skiwi_bool, skiwi_int64, skiwi_char_pointer, "(vertices->csv id \"file.csv\") exports the vertices of the object with tag `id` to a csv file.");
  register_external_primitive("view-bg-set!", (void*)&scm_set_bg_color, skiwi_void, skiwi_int64, skiwi_int64, skiwi_int64, "(view-bg-set! r g b) changes the background color to (r g b).");
  register_external_primitive("view-cs", (void*)&scm_get_view_coordinate_system, skiwi_scm, "(view-cs) returns the coordinate system of the view camera.");
  register_external_primitive("view-cs-set!", (void*)&scm_set_view_coordinate_system, skiwi_void, skiwi_scm, "(view-cs-set! cs) sets the coordinate system of the view camera. The coordinate system `cs` can be given as a vector of size 16 in column major format or as a list of lists in row major format.");
  register_external_primitive("view-edges-set!", (void*)&scm_set_edges, skiwi_void, skiwi_bool, "(view-edges-set! #t/#f) turns on/off rendering of edges.");
  register_external_primitive("view-export", (void*)&scm_export_image, skiwi_void, skiwi_char_pointer, "(view-export \"image-file.png\") exports the current view to a png image.");

  register_external_primitive("view-ear-left-detector-set!", (void*)&scm_set_show_ear_left_detector, skiwi_void, skiwi_bool, "(view-ear-left-detector-set! #t/#f) turns on/off rendering of the left ear detector result.");
  register_external_primitive("view-ear-right-detector-set!", (void*)&scm_set_show_ear_right_detector, skiwi_void, skiwi_bool, "(view-ear-right-detector-set! #t/#f) turns on/off rendering of the right ear detector result.");
  register_external_primitive("view-face-detector-set!", (void*)&scm_set_show_face_detector, skiwi_void, skiwi_bool, "(view-face-detector-set! #t/#f) turns on/off rendering of the face detector result.");

  register_external_primitive("view-hide!", (void*)&scm_hide_view, skiwi_void, "(view-hide!) hides the 3d view.");
  register_external_primitive("view-onebit-set!", (void*)&scm_set_one_bit, skiwi_void, skiwi_bool, "(view-onebit-set! #t/#f) turns on/off one-bit rendering.");
  register_external_primitive("view-position", (void*)&scm_get_position, skiwi_scm, skiwi_scm, skiwi_scm, "(view-position x y) returns the 3d position of the vertex in the last render of the 3d view at coordinate (x,y) .");
  register_external_primitive("view-index", (void*)&scm_get_index, skiwi_scm, skiwi_scm, skiwi_scm, "(view-index x y) returns the vertex index of the vertex in the last render of the 3d view at coordinate (x,y).");
  register_external_primitive("view-id", (void*)&scm_get_id, skiwi_scm, skiwi_scm, skiwi_scm, "(view-id x y) returns the id of the object in the last render of the 3d view at coordinate (x,y).");

  register_external_primitive("view-shading-set!", (void*)&scm_set_shading, skiwi_void, skiwi_bool, "(view-shading-set! #t/#f) turns on/off lighting.");
  register_external_primitive("view-shadow-set!", (void*)&scm_set_shadow, skiwi_void, skiwi_bool, "(view-shadow-set! #t/#f) turns on/off rendering of shadow.");

  register_external_primitive("view-show!", (void*)&scm_show_view, skiwi_void, "(view-show!) shows the 3d view.");
  register_external_primitive("view-size-set!", (void*)&scm_set_image_size, skiwi_void, skiwi_scm, skiwi_scm, "(view-size-set! w h) resizes the plotted image to size (w, h).");
  register_external_primitive("view-textured-set!", (void*)&scm_set_textured, skiwi_void, skiwi_bool, "(view-textured-set! #t/#f) turns on/off rendering of texture.");
  register_external_primitive("view-unzoom!", (void*)&scm_unzoom, skiwi_void, "(view-unzoom!) sets the camera to its initial position.");
  register_external_primitive("view-wireframe-set!", (void*)&scm_set_wireframe, skiwi_void, skiwi_bool, "(view-wireframe-set! #t/#f) turns on/off rendering of wireframe.");  

  register_external_primitive("exit", (void*)&scm_exit, skiwi_void, "(exit) can be used in the input script to end meshscript, so the REPL is skipped.");
  return nullptr;
  }

struct scheme_loop_data
  {
  scheme_loop_data() : initialised(false) {}
  std::mutex mt;
  std::condition_variable cv;
  bool initialised;
  std::unique_ptr<std::thread> t;
  };

std::string get_help_text()
  {
  std::string help = R"(This is meshscript. You are interacting with the REPL.
Enter scheme commands or one of the following:

,asm
,env
,exit
,expand
,external <optional: (part of) external function name>
,mem
,unresolved

)";
  return help;
  }

void create_scheme_with_loop(scheme_loop_data* sld, int argc, char** argv)
  {
  sld->mt.lock();
  sld->initialised = true;
  sld->cv.notify_all();
  sld->mt.unlock();

  skiwi::skiwi_parameters pars;
  pars.heap_size = 64 * 1024 * 1024;
  skiwi::set_prompt("ms> ");
  skiwi::set_welcome_message("\nWelcome to meshscript\nType ,? for help.\n");
  skiwi::set_help_text(get_help_text());
  skiwi::scheme_with_skiwi(&register_functions, nullptr, pars);

  for (int i = 1; i < argc; ++i)
    {
    std::string filename(argv[i]);
    skiwi::skiwi_runf(filename);
    }

  if (!quit)
    skiwi::skiwi_repl();

  g_view->quit();

  skiwi::skiwi_quit();
  }

std::unique_ptr<std::thread> create_threaded_scheme_loop(scheme_loop_data& sld, int argc, char** argv)
  {
  assert(!sld.initialised);
  std::unique_ptr<std::thread> res(new std::thread(create_scheme_with_loop, &sld, argc, argv));
  std::unique_lock<std::mutex> lk(sld.mt);
  if (!sld.initialised)
    sld.cv.wait(lk, [&] {return sld.initialised; });
  return res;
  }

void create_scheme_loop(scheme_loop_data& sld, int argc, char** argv)
  {
  sld.t = create_threaded_scheme_loop(sld, argc, argv);
  }

void close_scheme_loop(scheme_loop_data& sld)
  {
  sld.t->join();
  }

int main(int argc, char** argv)
  {
  if (SDL_Init(SDL_INIT_VIDEO))
    return -1;

  {
  view v;
  g_view = &v;
  scheme_loop_data sld;
  create_scheme_loop(sld, argc, argv);
  v.loop();
  close_scheme_loop(sld);
  }
  
  SDL_Quit();
  return 0;
  }