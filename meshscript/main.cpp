#ifndef _SDL_main_h
#define _SDL_main_h
#endif
#include <SDL.h>

#include <stdio.h>
#include <string.h>

#include <libskiwi/libskiwi.h>

#include <libskiwi/runtime.h>

#include "view.h"
#include "jet.h"
#include <thread>
#include <memory>
#include <mutex>
#include <sstream>
#include <iostream>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

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
      std::cout << "error: set-cs: second input parameter should be a vector of size 16\n";
      std::cout << "               current size of vector is " << v.size() << "\n";
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
    std::cout << "error: set-cs: invalid input type\n";
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

void set_color(int64_t id, int64_t clr_id)
  {
  g_view.v->set_color(id, clr_id);
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
  register_external_primitive("exit", &scm_exit, skiwi_void, "(exit) can be used in the input script to end meshscript");
  register_external_primitive("hide-view", &scm_hide_view, skiwi_void, "(hide-view) hides the 3d view");
  register_external_primitive("show-view", &scm_show_view, skiwi_void, "(show-view) shows the 3d view");  
  register_external_primitive("load-mesh", &load_mesh, skiwi_int64, skiwi_char_pointer, "(load-mesh \"stlfile.stl\") loads the stl file and returns an id. Similarly (load-mesh \"objfile.obj\") loads an obj file and return the id.");  
  register_external_primitive("load-pc", &load_pc, skiwi_int64, skiwi_char_pointer, "(load-pc \"pointcloud.ply\") loads the ply file as point cloud and returns an id.");
  register_external_primitive("set-cs", &scm_set_coordinate_system, skiwi_void, skiwi_int64, skiwi_scm, "(set-cs id cs) sets a new coordinate system for mesh id. The coordinate system cs can be given as a vector of size 16 in column major format or as a list of lists in row major format.");
  register_external_primitive("get-cs", &scm_get_coordinate_system, skiwi_scm, skiwi_int64, "(get-cs id) returns the coordinate system for mesh id.");
  register_external_primitive("make-mesh", &make_mesh, skiwi_int64, skiwi_scm, skiwi_scm, "(make-mesh vertices triangles) plots the mesh with given vertices and triangles, and returns the id of the plotted object. vertices should be a list of lists of the form ((x y z) (x y z) ...) with x,y,z floating point values, and triangles should be a list of list of the form ((a b c) (d e f) ...) with a,b... fixnums referring to the vertex indices.");
  register_external_primitive("vertices-to-csv", &vertices_to_csv, skiwi_bool, skiwi_int64, skiwi_char_pointer, "(vertices-to-csv id \"file.csv\") exports the vertices of mesh id to a csv file");
  register_external_primitive("triangles-to-csv", &triangles_to_csv, skiwi_bool, skiwi_int64, skiwi_char_pointer, "(triangles-to-csv id \"file.csv\") exports the triangles of mesh id to a csv file");
  register_external_primitive("show", &show, skiwi_void, skiwi_int64, "(show id) shows mesh id");
  register_external_primitive("hide", &hide, skiwi_void, skiwi_int64, "(hide id) hides mesh id");
  register_external_primitive("set-color", &set_color, skiwi_void, skiwi_int64, skiwi_int64, "(set-color id matcap-id) changes the matcap of mesh id. The matcap is given by its id matcap-id.");
  register_external_primitive("jet", &scm_jet, skiwi_scm, skiwi_scm, "(jet lst) takes a list of values between 0 and 1 and returns a list of lists with (r g b) values");
  register_external_primitive("set-vertex-colors", &scm_set_vertex_colors, skiwi_void, skiwi_int64, skiwi_scm, "(set-vertex-colors id clrlst) sets vertex colors for mesh id. The vertex colors are given as a list of lists with (r g b) values.");
  register_external_primitive("rotate", &scm_rotate, skiwi_void, skiwi_int64, skiwi_scm, skiwi_scm, skiwi_scm, "(rotate id x y z) rotates mesh id by x degrees over the x-axis, by y degrees over the y-axis, and by z degrees over the z_axis");
  register_external_primitive("translate", &scm_translate, skiwi_void, skiwi_int64, skiwi_scm, skiwi_scm, skiwi_scm, "(translate id x y z) translates mesh id by vector (x y z)");
  register_external_primitive("set-edges", &scm_set_edges, skiwi_void, skiwi_bool, "(set-edges #t/#f) turns on/off rendering of edges");
  register_external_primitive("set-shading", &scm_set_shading, skiwi_void, skiwi_bool, "(set-shading #t/#f) turns on/off lighting");
  register_external_primitive("set-textured", &scm_set_textured, skiwi_void, skiwi_bool, "(set-textured #t/#f) turns on/off rendering of texture");
  register_external_primitive("set-shadow", &scm_set_shadow, skiwi_void, skiwi_bool, "(set-shadow #t/#f) turns on/off rendering of shadow");
  register_external_primitive("set-wireframe", &scm_set_wireframe, skiwi_void, skiwi_bool, "(set-wireframe #t/#f) turns on/off rendering of wireframe");
  register_external_primitive("set-one-bit", &scm_set_one_bit, skiwi_void, skiwi_bool, "(set-one-bit #t/#f) turns on/off one-bit rendering");
  register_external_primitive("set-image-size", &scm_set_image_size, skiwi_void, skiwi_scm, skiwi_scm, "(set-image-size w h) resizes the plotted image to size (w, h)");
  register_external_primitive("unzoom", &scm_unzoom, skiwi_void, "(unzoom) sets the camera to its initial position");
  register_external_primitive("export-image", &scm_export_image, skiwi_void, skiwi_char_pointer, "(export-image \"image-file.png\") exports the current view to a png image");
  register_external_primitive("set-bg-color", &scm_set_bg_color, skiwi_void, skiwi_int64, skiwi_int64, skiwi_int64, "(set-bg-color r g b) changes the background color to (r g b).");
  register_external_primitive("get-position", &scm_get_position, skiwi_scm, skiwi_scm, skiwi_scm, "(get-position x y) returns the 3D position of coordinate (x,y).");
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