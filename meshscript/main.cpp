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
#include <sstream>
#include <iostream>

#include <jtk/geometry.h>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#include <tbb/enumerable_thread_specific.h>

using namespace jtk;

constexpr double pi = 3.141592653589793238462643383279;

struct view_data
  {
  view_data() : v(nullptr), initialised(false) {}
  view* v;
  std::mutex mt;
  std::condition_variable cv;
  bool initialised;
  std::unique_ptr<std::thread> t;
  };

void create_view_with_loop(view_data* vd)
  {
  vd->mt.lock();
  vd->v = new view();
  vd->initialised = true;
  vd->cv.notify_all();
  vd->mt.unlock();

  vd->v->loop();
  }

std::unique_ptr<std::thread> create_threaded_view(view_data& vd)
  {
  assert(!vd.initialised);
  assert(vd.v == nullptr);
  std::unique_ptr<std::thread> res(new std::thread(create_view_with_loop, &vd));
  std::unique_lock<std::mutex> lk(vd.mt);
  if (!vd.initialised)
    vd.cv.wait(lk, [&] {return vd.initialised; });
  return res;
  }

void create_view(view_data& vd)
  {
  vd.t = create_threaded_view(vd);
  }

void close_view(view_data& vd)
  {
  if (!vd.v)
    return;
  vd.v->quit();
  vd.t->join();
  delete vd.v;
  vd.v = nullptr;
  vd.initialised = false;
  }

view_data g_view;

bool quit = false;

void scm_exit()
  {
  quit = true;
  }

void scm_hide_view()
  {
  g_view.v->hide();
  }

void scm_show_view()
  {
  g_view.v->show();
  }

int64_t load_mesh(const char* filename)
  {
  int64_t id = g_view.v->load_mesh_from_file(filename);
  return id;
  }

int64_t load_pc(const char* filename)
  {
  int64_t id = g_view.v->load_pc_from_file(filename);
  return id;
  }

/*
Input can be a vector of size 16 in column major format,
or a list of lists in row major format like ((1 0 0 13) (0 1 0 12) (0 0 1 15) (0 0 0 1))
which can be read from a csv file.
*/
void scm_set_coordinate_system(int64_t id, skiwi::scm_type scheme_variable)
  {
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
      g_view.v->set_coordinate_system((uint32_t)id, m);
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
      g_view.v->set_coordinate_system((uint32_t)id, m);
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
void scm_set_view_coordinate_system(skiwi::scm_type scheme_variable)
  {
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
      g_view.v->set_coordinate_system(m);
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
      g_view.v->set_coordinate_system(m);
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
  auto cs = g_view.v->get_coordinate_system((uint32_t)id);
  std::vector<scm_type> lst;
  for (int r = 0; r < 4; ++r)
    {
    std::vector<scm_type> row;
    for (int c = 0; c < 4; ++c)
      {
      row.push_back(make_flonum(cs[r + 4*c]));
      }
    lst.push_back(make_list(row));
    }
  return make_list(lst);
  }

uint64_t scm_get_view_coordinate_system()
  {
  using namespace skiwi;
  auto cs = g_view.v->get_coordinate_system();
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
  tbb::enumerable_thread_specific< void* > thread_specific_scheme_context(nullptr);

  skiwi::skiwi_compiled_function_ptr marching_cubes_fun = nullptr;

  double marching_cubes_distance_fun(double x, double y, double z)
    {
    void*& local_context = thread_specific_scheme_context.local();
    if (!local_context)
      local_context = skiwi::skiwi_clone_context(skiwi::skiwi_get_context());
    skiwi::scm_type res = skiwi::skiwi_run_raw(marching_cubes_fun, local_context, x, y, z);
    return res.get_number();
    }

  }

int64_t scm_marching_cubes(skiwi::scm_type bb, skiwi::scm_type dim, skiwi::scm_type iso, skiwi::scm_type fun)
  {
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
    
      int64_t id = g_view.v->marching_cubes(bounding, width, height, depth, isovalue, &marching_cubes_distance_fun);

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

int64_t make_mesh(skiwi::scm_type scm_vertices, skiwi::scm_type scm_triangles)
  {
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
    int64_t id = g_view.v->load_mesh(vertices, triangles);
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
  bool res = g_view.v->vertices_to_csv(id, filename);
  return res;
  }

bool triangles_to_csv(int64_t id, const char* filename)
  {
  bool res = g_view.v->triangles_to_csv(id, filename);
  return res;
  }

void show(int64_t id)
  {
  g_view.v->show(id);
  }

void hide(int64_t id)
  {
  g_view.v->hide(id);
  }

void set_matcap(int64_t id, int64_t clr_id)
  {
  g_view.v->set_matcap(id, clr_id);
  }

void scm_unzoom()
  {
  g_view.v->unzoom();
  }

uint64_t scm_jet(skiwi::scm_type mask)
  {
  using namespace skiwi;
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

void scm_set_vertex_colors(int64_t id, skiwi::scm_type scm_colors)
  {
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
  g_view.v->set_vertex_colors((uint32_t)id, colors);
  }
  catch (std::runtime_error e)
    {
    std::cout << "error: set-vertex-colors: " << e.what() << "\n";
    }
  }

void scm_rotate(int64_t id, skiwi::scm_type x_axis, skiwi::scm_type y_axis, skiwi::scm_type z_axis)
  {
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
  g_view.v->premultiply_coordinate_system((uint32_t)id, rot);
  }

void scm_translate(int64_t id, skiwi::scm_type x_axis, skiwi::scm_type y_axis, skiwi::scm_type z_axis)
  {
  double x = x_axis.get_number();
  double y = y_axis.get_number();
  double z = z_axis.get_number();
  
  auto t = jtk::make_translation((float)x, (float)y, (float)z);
  g_view.v->premultiply_coordinate_system((uint32_t)id, t);
  }

void scm_set_shading(bool b)
  {
  g_view.v->set_shading(b);
  }

void scm_set_edges(bool b)
  {
  g_view.v->set_edges(b);
  }

void scm_set_textured(bool b)
  {
  g_view.v->set_textured(b);
  }

void scm_set_shadow(bool b)
  {
  g_view.v->set_shadow(b);
  }

void scm_set_wireframe(bool b)
  {
  g_view.v->set_wireframe(b);
  }

void scm_set_one_bit(bool b)
  {
  g_view.v->set_one_bit(b);
  }

void scm_set_image_size(skiwi::scm_type w, skiwi::scm_type h)
  {
  int width = (int)w.get_number();
  int height = (int)h.get_number();
  g_view.v->set_image_size(width, height);
  }

void scm_export_image(const char* filename)
  {
  auto im = g_view.v->get_image();
  stbi_write_png(filename, im.width(), im.height(), 4, (void*)im.data(), im.stride() * 4);
  }

void scm_set_bg_color(int64_t r, int64_t g, int64_t b)
  {
  g_view.v->set_bg_color((uint8_t)r, (uint8_t)g, (uint8_t)b);
  }

uint64_t scm_get_position(skiwi::scm_type x, skiwi::scm_type y)
  {
  using namespace skiwi;
  auto pos = g_view.v->get_world_position((int)x.get_number(), (int)y.get_number());
  std::vector<scm_type> coord;
  coord.push_back(make_flonum(pos[0]));
  coord.push_back(make_flonum(pos[1]));
  coord.push_back(make_flonum(pos[2]));
  return make_list(coord);
  }

void* register_functions(void*)
  {
  using namespace skiwi;
  register_external_primitive("cs-ref", &scm_get_coordinate_system, skiwi_scm, skiwi_int64, "(cs-ref id) returns the coordinate system for the object with tag `id`.");
  register_external_primitive("cs-rotate!", &scm_rotate, skiwi_void, skiwi_int64, skiwi_scm, skiwi_scm, skiwi_scm, "(cs-rotate! id x y z) rotates the object with tag `id` by `x` degrees over the x-axis, by `y` degrees over the y-axis, and by `z` degrees over the z_axis.");
  register_external_primitive("cs-set!", &scm_set_coordinate_system, skiwi_void, skiwi_int64, skiwi_scm, "(cs-set! id cs) sets a new coordinate system for the object with tag `id`. The coordinate system `cs` can be given as a vector of size 16 in column major format or as a list of lists in row major format.");
  register_external_primitive("cs-translate!", &scm_translate, skiwi_void, skiwi_int64, skiwi_scm, skiwi_scm, skiwi_scm, "(cs-translate! id x y z) translates the object with tag `id` by vector (x y z).");
  register_external_primitive("hide!", &hide, skiwi_void, skiwi_int64, "(hide! id) makes the object with tag `id` invisible.");
  register_external_primitive("jet", &scm_jet, skiwi_scm, skiwi_scm, "(jet lst) takes a list of values between 0 and 1 and returns a list of lists with (r g b) values.");
  register_external_primitive("load-mesh", &load_mesh, skiwi_int64, skiwi_char_pointer, "(load-mesh \"stlfile.stl\") loads the stl file and returns an id. Similarly (load-mesh \"objfile.obj\") loads an obj file and returns the id.");
  register_external_primitive("load-pointcloud", &load_pc, skiwi_int64, skiwi_char_pointer, "(load-pointcloud \"pointcloud.ply\") loads the ply file as point cloud and returns an id.");
  register_external_primitive("make-mesh", &make_mesh, skiwi_int64, skiwi_scm, skiwi_scm, "(make-mesh vertices triangles) plots the mesh with given vertices and triangles, and returns the id of the plotted object. Vertices should be a list of lists of the form ((x y z) (x y z) ...) with x,y,z floating point values, and triangles should be a list of lists of the form ((a b c) (d e f) ...) with a,b... fixnums referring to the vertex indices.");
  register_external_primitive("marching-cubes", &scm_marching_cubes, skiwi_int64, skiwi_scm, skiwi_scm, skiwi_scm, skiwi_scm, "(marching-cubes bb dim isovalue fun) with bb of the form ((min_x max_x) (min_y max_y) (min_z max_z)), dim of the form (width height depth), isovalue a flonum, fun a lambda function accepting (x y z) values and returning a distance.");
  register_external_primitive("matcap-set!", &set_matcap, skiwi_void, skiwi_int64, skiwi_int64, "(matcap-set! id matcap-id) changes the matcap of the object with tag `id`. The matcap is given by its id matcap-id.");  
  register_external_primitive("show!", &show, skiwi_void, skiwi_int64, "(show! id) makes the object with tag `id` visible.");
  register_external_primitive("triangles->csv", &triangles_to_csv, skiwi_bool, skiwi_int64, skiwi_char_pointer, "(triangles->csv id \"file.csv\") exports the triangles of the object with tag `id` to a csv file.");
  register_external_primitive("vertexcolors-set!", &scm_set_vertex_colors, skiwi_void, skiwi_int64, skiwi_scm, "(vertexcolors-set! id clrlst) sets vertex colors for the object with tag `id`. The vertex colors are given as a list of lists with (r g b) values.");
  register_external_primitive("vertices->csv", &vertices_to_csv, skiwi_bool, skiwi_int64, skiwi_char_pointer, "(vertices->csv id \"file.csv\") exports the vertices of the object with tag `id` to a csv file.");
  register_external_primitive("view-bg-set!", &scm_set_bg_color, skiwi_void, skiwi_int64, skiwi_int64, skiwi_int64, "(view-bg-set! r g b) changes the background color to (r g b).");
  register_external_primitive("view-cs", &scm_get_view_coordinate_system, skiwi_scm, "(view-cs) returns the coordinate system of the view camera.");
  register_external_primitive("view-cs-set!", &scm_set_view_coordinate_system, skiwi_void, skiwi_scm, "(view-cs-set! cs) sets the coordinate system of the view camera.");
  register_external_primitive("view-edges-set!", &scm_set_edges, skiwi_void, skiwi_bool, "(view-edges-set! #t/#f) turns on/off rendering of edges.");
  register_external_primitive("view-export", &scm_export_image, skiwi_void, skiwi_char_pointer, "(view-export \"image-file.png\") exports the current view to a png image.");
  register_external_primitive("view-hide!", &scm_hide_view, skiwi_void, "(view-hide!) hides the 3d view.");
  register_external_primitive("view-onebit-set!", &scm_set_one_bit, skiwi_void, skiwi_bool, "(view-onebit-set! #t/#f) turns on/off one-bit rendering.");
  register_external_primitive("view-ref", &scm_get_position, skiwi_scm, skiwi_scm, skiwi_scm, "(view-ref x y) returns the 3D position of coordinate (x,y).");
  register_external_primitive("view-shading-set!", &scm_set_shading, skiwi_void, skiwi_bool, "(view-shading-set! #t/#f) turns on/off lighting.");
  register_external_primitive("view-shadow-set!", &scm_set_shadow, skiwi_void, skiwi_bool, "(view-shadow-set! #t/#f) turns on/off rendering of shadow.");
  register_external_primitive("view-show!", &scm_show_view, skiwi_void, "(view-show!) shows the 3d view.");  
  register_external_primitive("view-size-set!", &scm_set_image_size, skiwi_void, skiwi_scm, skiwi_scm, "(view-size-set! w h) resizes the plotted image to size (w, h).");  
  register_external_primitive("view-textured-set!", &scm_set_textured, skiwi_void, skiwi_bool, "(view-textured-set! #t/#f) turns on/off rendering of texture.");
  register_external_primitive("view-unzoom!", &scm_unzoom, skiwi_void, "(view-unzoom!) sets the camera to its initial position.");
  register_external_primitive("view-wireframe-set!", &scm_set_wireframe, skiwi_void, skiwi_bool, "(view-wireframe-set! #t/#f) turns on/off rendering of wireframe.");
 
  register_external_primitive("exit", &scm_exit, skiwi_void, "(exit) can be used in the input script to end meshscript.");
  return nullptr;
  }

int main(int argc, char** argv)
  {
  if (SDL_Init(SDL_INIT_VIDEO))
    return -1;

  create_view(g_view);

  skiwi::skiwi_parameters pars;
  pars.heap_size = 64 * 1024 * 1024;
  skiwi::set_prompt("ms> ");
  skiwi::scheme_with_skiwi(&register_functions, nullptr, pars);

  for (int i = 1; i < argc; ++i)
    {
    std::string filename(argv[i]);
    skiwi::skiwi_runf(filename);
    }

  if (!quit)
    skiwi::skiwi_repl();

  skiwi::skiwi_quit();

  close_view(g_view);
  SDL_Quit();
  return 0;
  }