#include "view.h"
#include "mesh.h"
#include "mm.h"
#include "pc.h"
#include "sp.h"
#include "im.h"
#include "ear_detector.h"
#include "face_detector.h"
#include "fill_holes.h"
#include "view.h"
#include "distance_map.h"
#include "shape_predictor.h"
#include "lscm.h"

#include <libpoisson/poisson_reconstruction_screened.h>

#include <libcork/cork.h>

#include <iostream>

#include <SDL_syswm.h>

#include <jtk/geometry.h>

#include <icp/icp_point_to_plane.h>

#include <sstream>

using namespace jtk;

namespace
  {

  uint32_t get_next_power_of_two(uint32_t v)
    {
    --v;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    ++v;
    return v;
    }

  void gl_check_error(const char* txt)
    {
    unsigned int err = glGetError();
    if (err)
      {
      std::stringstream str;
      str << "GL error " << err << ": " << txt;
      throw std::runtime_error(str.str());
      }
    }

  void clear_screen(image<uint32_t>& screen)
    {
    for (auto& v : screen)
      v = 0xff000000 | (uint32_t(49) << 16) | (uint32_t(49) << 8) | uint32_t(49);
    }

  }

view::view() : _w(1600), _h(900), _window(nullptr)
  {
  SDL_DisplayMode DM;
  SDL_GetCurrentDisplayMode(0, &DM);
  _w_max = DM.w;
  _h_max = DM.h;

  if (_w > _w_max)
    _w = _w_max;
  if (_h > _h_max)
    _h = _h_max;
  // Setup window
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 8);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
  SDL_DisplayMode current;
  SDL_GetCurrentDisplayMode(0, &current);

  _quit = false;
  _refresh = true;
  _show_ear_left_detector = false;
  _show_ear_right_detector = false;
  _show_face_detector = false;

  //prepare_window();

  _screen = image<uint32_t>(_w, _h);
  clear_screen(_screen);

  _m.left_dragging = false;
  _m.right_dragging = false;
  _m.right_button_down = false;
  _m.left_button_down = false;
  _m.wheel_down = false;
  _m.wheel_mouse_pressed = false;
  _m.ctrl_pressed = false;
  _m.mouse_x = 0.f;
  _m.mouse_y = 0.f;
  _m.prev_mouse_x = 0.f;
  _m.prev_mouse_y = 0.f;
  _m.wheel_rotation = 0.f;

  _canvas.resize(800, 600);
  _canvas_pos_x = ((int32_t)_w - (int32_t)_canvas.width()) / 2;
  if (_canvas_pos_x & 3)
    _canvas_pos_x += 4 - (_canvas_pos_x & 3);
  _canvas_pos_y = ((int32_t)_h - (int32_t)_canvas.height()) / 2;

  _settings.one_bit = false;
  _settings.shadow = false;
  _settings.edges = true;
  _settings.wireframe = false;
  _settings.shading = true;
  _settings.textured = true;
  _settings.vertexcolors = true;

  _suspend = false;
  _resume = false;

  p_face_detector.reset(new face_detector());
  p_ear_detector.reset(new ear_detector());
  }

view::~view()
  {
  delete_window();
  }

void view::delete_window()
  {
  if (!_window)
    return;
  glDeleteTextures(1, &_gl_texture);
  SDL_DestroyWindow(_window);
  _window = nullptr;
  }

void view::prepare_window()
  {
  if (_window)
    return;
  _window = SDL_CreateWindow("meshscript", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, _w, _h, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
      
  if (_window == NULL)
  {
      std::cout << SDL_GetError() << "\n";
      return;
  }
  SDL_GLContext gl_context = SDL_GL_CreateContext(_window);
  SDL_GL_SetSwapInterval(1); // Enable vsync

  SDL_GL_MakeCurrent(_window, gl_context);

  glEnable(GL_TEXTURE_2D);
  glGenTextures(1, &_gl_texture);
  glBindTexture(GL_TEXTURE_2D, _gl_texture);

  _gl_texture_w = get_next_power_of_two(_w_max);
  _gl_texture_h = get_next_power_of_two(_h_max);

  GLint max_texture_size;
  glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max_texture_size);

  if (_gl_texture_w > (uint32_t)max_texture_size)
    _gl_texture_w = (uint32_t)max_texture_size;
  if (_gl_texture_h > (uint32_t)max_texture_size)
    _gl_texture_h = (uint32_t)max_texture_size;

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, _gl_texture_w, _gl_texture_h, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
  gl_check_error("glTexImage2D in view.cpp");
  }

int64_t view::make_cube(float w, float h, float d)
  {
  std::scoped_lock lock(_mut);
  mesh* new_object;
  uint32_t id;
  _db.create_mesh(new_object, id);
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  ::make_cube(*new_object, w, h, d);
  if (new_object->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return id;
  }

int64_t view::make_sphere(float r, uint32_t subdivision_levels)
  {
  std::scoped_lock lock(_mut);
  mesh* new_object;
  uint32_t id;
  _db.create_mesh(new_object, id);
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  ::make_sphere(*new_object, r, subdivision_levels);
  if (new_object->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return id;
  }

int64_t view::make_icosahedron(float r)
  {
  std::scoped_lock lock(_mut);
  mesh* new_object;
  uint32_t id;
  _db.create_mesh(new_object, id);
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  ::make_icosahedron(*new_object, r);
  if (new_object->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return id;
  }

int64_t view::make_cylinder(float r, float d, uint32_t n)
  {
  std::scoped_lock lock(_mut);
  mesh* new_object;
  uint32_t id;
  _db.create_mesh(new_object, id);
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  ::make_cylinder(*new_object, r, d, n);
  if (new_object->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return id;
  }

int64_t view::triangulate(const std::vector<jtk::vec2<float>>& vertices)
  {
  std::scoped_lock lock(_mut);
  mesh* new_object;
  uint32_t id;
  _db.create_mesh(new_object, id);
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  ::triangulate_points(*new_object, vertices);
  if (new_object->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return id;
  }

int64_t view::extrude(const std::vector<jtk::vec2<float>>& vertices, float h)
  {
  std::scoped_lock lock(_mut);
  mesh* new_object;
  uint32_t id;
  _db.create_mesh(new_object, id);
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  ::extrude_points(*new_object, vertices, h);
  if (new_object->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return id;
  }

int64_t view::revolve(const std::vector<jtk::vec2<float>>& vertices, uint32_t n, bool closed)
  {
  std::scoped_lock lock(_mut);
  mesh* new_object;
  uint32_t id;
  _db.create_mesh(new_object, id);
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  ::revolve_points(*new_object, vertices, n, closed);
  if (new_object->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return id;
  }

void view::subdivide(uint32_t id, uint32_t nr)
  {
  std::scoped_lock lock(_mut);
  std::vector<jtk::vec3<float>>* p_vertices = get_vertices(_db, id);
  std::vector<jtk::vec3<uint32_t>>* p_triangles = get_triangles(_db, id);
  if (p_vertices && p_triangles)
    {
    for (uint32_t i = 0; i < nr; ++i)
      jtk::dyadic_subdivide(*p_vertices, *p_triangles);
    remove_object(id, _scene);
    if (is_visible(_db, id))
      add_object(id, _scene, _db);
    _refresh = true;
    }
  }

void view::smooth(uint32_t id, uint32_t iterations, float lambda, float mu)
  {
  std::scoped_lock lock(_mut);
  std::vector<jtk::vec3<float>>* p_vertices = get_vertices(_db, id);
  std::vector<jtk::vec3<uint32_t>>* p_triangles = get_triangles(_db, id);
  if (p_vertices && p_triangles)
    {
    jtk::smooth(*p_vertices, *p_triangles, iterations, lambda, mu);
    remove_object(id, _scene);
    if (is_visible(_db, id))
      add_object(id, _scene, _db);
    _refresh = true;
    }
  }

void view::scale(uint32_t id, float sx, float sy, float sz)
  {
  std::scoped_lock lock(_mut);
  std::vector<jtk::vec3<float>>* p_vertices = get_vertices(_db, id);
  if (p_vertices)
    {
    for (auto& v : *p_vertices)
      {
      v[0] *= sx;
      v[1] *= sy;
      v[2] *= sz;
      }
    remove_object(id, _scene);
    if (is_visible(_db, id))
      add_object(id, _scene, _db);
    prepare_scene(_scene);
    _refresh = true;
    }
  }

int64_t view::duplicate(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh(id);
  pc* p = _db.get_pc(id);
  mm* m2 = _db.get_mm(id);
  sp* s = _db.get_sp(id);
  im* i = _db.get_image(id);
  if (m)
    {
    mesh* new_object;
    uint32_t id_out;
    _db.create_mesh(new_object, id_out);
    _matcap.map_db_id_to_matcap_id(id_out, _get_semirandom_matcap_id(id_out));
    *new_object = *m;
    if (new_object->visible)
      add_object(id_out, _scene, _db);
    prepare_scene(_scene);    
    _refresh = true;
    return id_out;
    }
  if (p)
    {
    pc* new_object;
    uint32_t id_out;
    _db.create_pc(new_object, id_out);
    *new_object = *p;
    if (new_object->visible)
      add_object(id_out, _scene, _db);
    prepare_scene(_scene);
    _refresh = true;
    return id_out;
    }
  if (m2)
    {
    mm* new_object;
    uint32_t id_out;
    _db.create_mm(new_object, id_out);
    _matcap.map_db_id_to_matcap_id(id_out, _get_semirandom_matcap_id(id_out));
    *new_object = *m2;
    if (new_object->visible)
      add_object(id_out, _scene, _db);
    prepare_scene(_scene);
    _refresh = true;
    return id_out;
    }
  if (s)
    {
    sp* new_object;
    uint32_t id_out;
    _db.create_sp(new_object, id_out);
    *new_object = *s;   
    return id_out;
    }
  if (i)
    {
    im* new_object;
    uint32_t id_out;
    _db.create_image(new_object, id_out);
    *new_object = *i;
    return id_out;
    }
  return -1;
  }

int64_t view::mesh_texture_to_image(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh(id);
  if (m)
    {
    im* db_image;
    uint32_t id;
    _db.create_image(db_image, id);
    db_image->texture = m->texture;
    if (db_image->texture.width() == 0 || db_image->texture.height() == 0)
      {
      db_image->texture = make_dummy_texture(512, 512);
      }
    return id;
    }
  return -1;
  }

int64_t view::mesh_to_pointcloud(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh(id);
  if (m)
    {
    pc* db_pc;
    uint32_t id_out;
    _db.create_pc(db_pc, id_out);

    std::vector<vec3<float>> triangle_normals, vertex_normals;
    compute_triangle_normals(triangle_normals, m->vertices.data(), m->triangles.data(), m->triangles.size());
    compute_vertex_normals(vertex_normals, triangle_normals.data(), m->vertices.data(), m->vertices.size(), m->triangles.data(), m->triangles.size());
    db_pc->vertices = m->vertices;       
    db_pc->vertex_colors = convert_vertex_colors(m->vertex_colors);   
    db_pc->normals = vertex_normals;
    db_pc->cs = m->cs;
    db_pc->visible = true;
    if (db_pc->visible)
      add_object(id_out, _scene, _db);
    prepare_scene(_scene);
    ::unzoom(_scene);
    _refresh = true;
    return id_out;
    }
  return -1;
  }

int64_t view::load_mesh_from_file(const char* filename)
  {
  std::scoped_lock lock(_mut);
  mesh m;
  std::string f(filename);
  bool res = read_from_file(m, f);
  if (!res)
    return -1;
  mesh* db_mesh;
  uint32_t id;
  _db.create_mesh(db_mesh, id);
  db_mesh->vertices.swap(m.vertices);
  db_mesh->triangles.swap(m.triangles);
  db_mesh->uv_coordinates.swap(m.uv_coordinates);
  db_mesh->texture.swap(m.texture);
  db_mesh->vertex_colors.swap(m.vertex_colors);
  db_mesh->cs = m.cs;
  db_mesh->visible = m.visible;
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  if (db_mesh->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return (int64_t)id;
  }

int64_t view::load_morphable_model_from_file(const char* filename)
  {
  std::scoped_lock lock(_mut);
  mm m;
  std::string f(filename);
  bool res = read_from_file(m, f);
  if (!res)
    return -1;
  mm* db_mm;
  uint32_t id;
  _db.create_mm(db_mm, id);
  swap(db_mm->shape, m.shape);
  swap(db_mm->color, m.color);
  db_mm->vertices.swap(m.vertices);
  db_mm->coefficients.swap(m.coefficients);
  db_mm->color_coefficients.swap(m.color_coefficients);
  db_mm->vertex_colors.swap(m.vertex_colors);
  db_mm->cs = m.cs;
  db_mm->visible = m.visible;
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  if (db_mm->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return (int64_t)id;
  }

int64_t view::load_pc_from_file(const char* filename)
  {
  std::scoped_lock lock(_mut);
  pc point_cloud;
  std::string f(filename);
  bool res = read_from_file(point_cloud, f);
  if (!res)
    return -1;
  pc* db_pc;
  uint32_t id;
  _db.create_pc(db_pc, id);
  db_pc->vertices.swap(point_cloud.vertices);
  db_pc->normals.swap(point_cloud.normals);
  db_pc->vertex_colors.swap(point_cloud.vertex_colors);
  db_pc->cs = point_cloud.cs;
  db_pc->visible = point_cloud.visible;
  if (db_pc->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return (int64_t)id;
  }

int64_t view::load_image_from_file(const char* filename)
  {
  std::scoped_lock lock(_mut);
  im image;
  std::string f(filename);
  bool res = read_from_file(image, f);
  if (!res)
    return -1;
  im* db_image;
  uint32_t id;
  _db.create_image(db_image, id);
  db_image->texture.swap(image.texture);  
  return (int64_t)id;
  }

int64_t view::load_mesh(const std::vector<vec3<float>>& vertices, const std::vector<vec3<uint32_t>>& triangles)
  {
  std::scoped_lock lock(_mut);
  mesh* db_mesh;
  uint32_t id;
  _db.create_mesh(db_mesh, id);
  db_mesh->vertices = vertices;
  db_mesh->triangles = triangles;
  db_mesh->cs = get_identity();
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  db_mesh->visible = true;
  add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return (int64_t)id;
  }

int64_t view::load_pointcloud(const std::vector<vec3<float>>& vertices)
  {
  std::scoped_lock lock(_mut);
  pc* db_pc;
  uint32_t id;
  _db.create_pc(db_pc, id);
  db_pc->vertices = vertices;
  db_pc->cs = get_identity();
  db_pc->visible = true;
  add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return (int64_t)id;
  }

void view::set_coordinate_system(uint32_t id, const float4x4& cs)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh(id);
  if (m)
    {
    m->cs = cs;
    remove_object(id, _scene);
    if (m->visible)
      add_object(id, _scene, _db);
    }
  pc* p = _db.get_pc(id);
  if (p)
    {
    p->cs = cs;
    remove_object(id, _scene);
    if (p->visible)
      add_object(id, _scene, _db);
    }
  mm* mo = _db.get_mm((uint32_t)id);
  if (mo)
    {
    mo->cs = cs;
    remove_object(id, _scene);
    if (mo->visible)
      add_object(id, _scene, _db);
    }
  if (m || p || mo)
    {
    prepare_scene(_scene);
    _refresh = true;
    }
  }


void view::set_coordinate_system(const float4x4& cs)
  {
  std::scoped_lock lock(_mut);
  _scene.coordinate_system = cs;
  _scene.coordinate_system_inv = jtk::invert_orthonormal(cs);
  _refresh = true;
  }

void view::render_scene()
  {
  // assumes a lock has been set already
  _canvas.update_settings(_settings);
  _canvas.render_scene(&_scene);
  _pixels = _canvas.get_pixels();
  _canvas.canvas_to_image(_pixels, _matcap);
  _canvas.render_pointclouds_on_image(&_scene, _pixels);
  _refresh = false;

  if (_show_face_detector && p_face_detector.get())
    {
    auto rectangles = p_face_detector->detect(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data());
    p_face_detector->draw_prediction_rgba(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), rectangles);
    for (const auto& sp : _db.get_sps())
      {
      if (sp.second->p_shape_predictor.get() && sp.second->odl == sp::odl_facial)
        {
        for (const auto& r : rectangles)
          {
          auto points = sp.second->p_shape_predictor->predict(r, _canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), sp.second->flip_horizontal);
          sp.second->p_shape_predictor->draw_prediction_rgba(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), points);
          }
        }
      }
    }

  if (_show_ear_right_detector && p_ear_detector.get())
    {
    auto rectangles = p_ear_detector->detect(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), true);
    p_ear_detector->draw_prediction_rgba(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), rectangles);
    for (const auto& sp : _db.get_sps())
      {
      if (sp.second->p_shape_predictor.get() && sp.second->odl == sp::odl_ear_right)
        {
        for (const auto& r : rectangles)
          {
          auto points = sp.second->p_shape_predictor->predict(r, _canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), sp.second->flip_horizontal);
          sp.second->p_shape_predictor->draw_prediction_rgba(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), points);
          }
        }
      }
    }

  if (_show_ear_left_detector && p_ear_detector.get())
    {
    auto rectangles = p_ear_detector->detect(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), false);
    p_ear_detector->draw_prediction_rgba(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), rectangles);
    for (const auto& sp : _db.get_sps())
      {
      if (sp.second->p_shape_predictor.get() && sp.second->odl == sp::odl_ear_left)
        {
        for (const auto& r : rectangles)
          {
          auto points = sp.second->p_shape_predictor->predict(r, _canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), sp.second->flip_horizontal);
          sp.second->p_shape_predictor->draw_prediction_rgba(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), points);
          }
        }
      }
    }

  }

jtk::image<uint32_t> view::get_image()
  {
  std::scoped_lock lock(_mut);
  if (_refresh)
    {
    render_scene();
    }
  return _canvas.get_image();
  }

jtk::float4x4 view::get_coordinate_system(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh(id);
  if (m)
    return m->cs;
  pc* p = _db.get_pc(id);
  if (p)
    return p->cs;
  mm* mo = _db.get_mm((uint32_t)id);
  if (mo)
    return mo->cs;
  return jtk::get_identity();
  }

jtk::float4x4 view::get_coordinate_system()
  {
  std::scoped_lock lock(_mut);
  return _scene.coordinate_system;
  }

jtk::float4x4 view::icp(uint32_t id1, uint32_t id2, double inlier_distance)
  {
  std::scoped_lock lock(_mut);
  std::vector<jtk::vec3<float>> model_points, template_points;
  auto vert_id1 = get_vertices(_db, id1);
  auto vert_id2 = get_vertices(_db, id2);
  if (!vert_id1 || !vert_id2)
    return jtk::get_identity();
  auto cs1 = *get_cs(_db, id1);
  auto cs2 = *get_cs(_db, id2);
  model_points = *vert_id1;
  template_points = *vert_id2;
  for (auto& v : model_points)
    {
    jtk::float4 V(v[0], v[1], v[2], 1.f);
    V = jtk::matrix_vector_multiply(cs1, V);
    v[0] = V[0];
    v[1] = V[1];
    v[2] = V[2];
    }
  for (auto& v : template_points)
    {
    jtk::float4 V(v[0], v[1], v[2], 1.f);
    V = jtk::matrix_vector_multiply(cs2, V);
    v[0] = V[0];
    v[1] = V[1];
    v[2] = V[2];
    }
  std::cout << "Initializing icp\n";
  jtk::icp_point_to_plane icp_algo(model_points);
  matf16 transf = jtk::identity<float>(4, 4);
  icp_algo.set_maximum_iterations(10);
  std::cout << "Running icp\n";
  float residu = icp_algo.fit(transf, template_points, (float)inlier_distance);
  std::cout << "icp residu: " << residu << "\n";
  jtk::float4x4 out;
  for (int i = 0; i < 16; ++i)
    {
    out[i] = transf(i % 4, i / 4);
    }
  return out;
  }

std::vector<float> view::distance_map(uint32_t id1, uint32_t id2, bool sign)
  {
  std::scoped_lock lock(_mut);
  std::vector<float> out;
  auto vert_id1 = get_vertices(_db, id1);
  auto vert_id2 = get_vertices(_db, id2);
  auto tria_id2 = get_triangles(_db, id2);
  if (!vert_id1)
    return out;
  if (!vert_id2)
    return out;
  if (!tria_id2)
    return out;
  auto cs1 = *get_cs(_db, id1);
  auto cs2 = *get_cs(_db, id2);
  out = ::distance_map(*vert_id1, cs1, *vert_id2, *tria_id2, cs2, sign);

  return out;
  }

void view::set_bg_color(uint8_t r, uint8_t g, uint8_t b)
  {
  std::scoped_lock lock(_mut);
  _canvas.set_background_color(r, g, b);
  _refresh = true;
  }

void view::premultiply_coordinate_system(uint32_t id, const jtk::float4x4& cs)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh(id);
  if (m)
    {
    m->cs = jtk::matrix_matrix_multiply(cs, m->cs);
    remove_object(id, _scene);
    if (m->visible)
      add_object(id, _scene, _db);
    }
  pc* p = _db.get_pc(id);
  if (p)
    {
    p->cs = jtk::matrix_matrix_multiply(cs, p->cs);
    remove_object(id, _scene);
    if (p->visible)
      add_object(id, _scene, _db);
    }
  mm* mo = _db.get_mm(id);
  if (mo)
    {
    mo->cs = jtk::matrix_matrix_multiply(cs, mo->cs);
    remove_object(id, _scene);
    if (mo->visible)
      add_object(id, _scene, _db);
    }
  if (m || p || mo)
    {
    prepare_scene(_scene);
    _refresh = true;
    }
  }

void view::mesh_texture_set(uint32_t id, uint32_t tex)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh(id);
  im* i = _db.get_image(tex);
  if (m && i)
    {
    m->texture = i->texture;
    remove_object(id, _scene);
    if (m->visible)
      add_object(id, _scene, _db);
    _refresh = true;
    }
  }

void view::set_vertex_colors(uint32_t id, const std::vector<jtk::vec3<uint8_t>>& colors)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh(id);
  if (m)
    {
    m->vertex_colors.clear();
    m->vertex_colors.reserve(colors.size());
    for (auto& clr : colors)
      m->vertex_colors.emplace_back(clr[0] / 255.f, clr[1] / 255.f, clr[2] / 255.f);
    remove_object(id, _scene);
    if (m->visible)
      add_object(id, _scene, _db);
    _refresh = true;
    }
  pc* p = _db.get_pc(id);
  if (p)
    {
    p->vertex_colors.clear();
    p->vertex_colors.reserve(colors.size());
    for (auto& clr : colors)
      p->vertex_colors.emplace_back(make_color(clr[0], clr[1], clr[2]));
    remove_object(id, _scene);
    if (p->visible)
      add_object(id, _scene, _db);
    _refresh = true;
    }
  }

bool view::vertices_to_csv(int64_t id, const char* filename)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (m)
    {
    std::string fn(filename);
    return ::vertices_to_csv(*m, fn);
    }
  pc* p = _db.get_pc((uint32_t)id);
  if (p)
    {
    std::string fn(filename);
    return ::vertices_to_csv(*p, fn);
    }
  mm* mo = _db.get_mm((uint32_t)id);
  if (mo)
    {
    std::string fn(filename);
    return ::vertices_to_csv(*mo, fn);
    }
  return false;
  }

bool view::triangles_to_csv(int64_t id, const char* filename)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (m)
    {
    std::string fn(filename);
    return ::triangles_to_csv(*m, fn);
    }
  mm* mo = _db.get_mm((uint32_t)id);
  if (mo)
    {
    std::string fn(filename);
    return ::triangles_to_csv(*mo, fn);
    }
  return false;
  }

void view::hide(int64_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  pc* p = _db.get_pc((uint32_t)id);
  mm* mo = _db.get_mm((uint32_t)id);
  if (m)
    m->visible = false;
  else if (p)
    p->visible = false;
  else if (mo)
    mo->visible = false;
  else
    return;
  remove_object((uint32_t)id, _scene);
  prepare_scene(_scene);
  _refresh = true;
  }

void view::hide()
  {
  std::scoped_lock lock(_mut);
  _suspend = true;
  }

void view::show()
  {
  std::scoped_lock lock(_mut);
  _resume = true;
  }

void view::set_matcap(uint32_t id, uint32_t clr_id)
  {
  std::scoped_lock lock(_mut);
  im* i = _db.get_image(clr_id);
  if (i)
    {
    matcap new_matcap;
    new_matcap.im = i->texture;
    new_matcap.cavity_clr = i->texture(0, 0);
    _matcap.make_new_matcap(clr_id, new_matcap);
    _matcap.map_db_id_to_matcap_id(id, clr_id);
    }
  else
    _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(clr_id));
  _refresh = true;
  }

void view::show(int64_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  pc* p = _db.get_pc((uint32_t)id);
  mm* mo = _db.get_mm((uint32_t)id);
  if (m)
    m->visible = true;
  else if (p)
    p->visible = true;
  else if (mo)
    mo->visible = true;
  else
    return;
  if (m || p || mo)
    {
    remove_object((uint32_t)id, _scene);
    add_object((uint32_t)id, _scene, _db);
    prepare_scene(_scene);
    _refresh = true;
    }
  }

void view::set_shading(bool b)
  {
  std::scoped_lock lock(_mut);
  _settings.shading = b;
  _refresh = true;
  }

void view::set_edges(bool b)
  {
  std::scoped_lock lock(_mut);
  _settings.edges = b;
  _refresh = true;
  }

void view::set_textured(bool b)
  {
  std::scoped_lock lock(_mut);
  _settings.textured = b;
  _refresh = true;
  }

void view::set_vertexcolors(bool b)
  {
  std::scoped_lock lock(_mut);
  _settings.vertexcolors = b;
  _refresh = true;
  }

void view::set_shadow(bool b)
  {
  std::scoped_lock lock(_mut);
  _settings.shadow = b;
  _refresh = true;
  }

void view::set_wireframe(bool b)
  {
  std::scoped_lock lock(_mut);
  _settings.wireframe = b;
  _refresh = true;
  }

void view::set_one_bit(bool b)
  {
  std::scoped_lock lock(_mut);
  _settings.one_bit = b;
  _refresh = true;
  }

void view::set_image_size(int w, int h)
  {
  if (w < 0 || h < 0)
    return;
  std::scoped_lock lock(_mut);
  resize_canvas(w, h);
  _refresh = true;
  }

void view::force_redraw()
  {
  std::scoped_lock lock(_mut);
  render_scene();
  }

jtk::vec3<float> view::get_world_position(int x, int y)
  {
  std::scoped_lock lock(_mut);
  jtk::vec3<float> invalid_vertex(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
  if (x < 0 || y < 0 || x >= (int)_pixels.width() || y >= (int)_pixels.height())
    return invalid_vertex;
  const auto& p = _pixels(x, y);
  if (p.db_id == 0)
    return invalid_vertex;
  mesh* m = _db.get_mesh((uint32_t)p.db_id);
  if (m)
    {
    const uint32_t v0 = m->triangles[p.object_id][0];
    const uint32_t v1 = m->triangles[p.object_id][1];
    const uint32_t v2 = m->triangles[p.object_id][2];
    const float4 V0(m->vertices[v0][0], m->vertices[v0][1], m->vertices[v0][2], 1.f);
    const float4 V1(m->vertices[v1][0], m->vertices[v1][1], m->vertices[v1][2], 1.f);
    const float4 V2(m->vertices[v2][0], m->vertices[v2][1], m->vertices[v2][2], 1.f);
    const float4 pos = V0 * (1.f - p.barycentric_u - p.barycentric_v) + p.barycentric_u*V1 + p.barycentric_v*V2;
    auto world_pos = matrix_vector_multiply(m->cs, pos);
    return jtk::vec3<float>(world_pos[0], world_pos[1], world_pos[2]);
    }
  pc* ptcl = _db.get_pc((uint32_t)p.db_id);
  if (ptcl)
    {
    const float4 pos(ptcl->vertices[p.object_id][0], ptcl->vertices[p.object_id][1], ptcl->vertices[p.object_id][2], 1.f);
    auto world_pos = matrix_vector_multiply(m->cs, pos);
    return jtk::vec3<float>(world_pos[0], world_pos[1], world_pos[2]);
    }
  mm* mo = _db.get_mm((uint32_t)p.db_id);
  if (mo)
    {
    const uint32_t v0 = mo->shape.triangles[p.object_id][0];
    const uint32_t v1 = mo->shape.triangles[p.object_id][1];
    const uint32_t v2 = mo->shape.triangles[p.object_id][2];
    const float4 V0(mo->vertices[v0][0], mo->vertices[v0][1], mo->vertices[v0][2], 1.f);
    const float4 V1(mo->vertices[v1][0], mo->vertices[v1][1], mo->vertices[v1][2], 1.f);
    const float4 V2(mo->vertices[v2][0], mo->vertices[v2][1], mo->vertices[v2][2], 1.f);
    const float4 pos = V0 * (1.f - p.barycentric_u - p.barycentric_v) + p.barycentric_u*V1 + p.barycentric_v*V2;
    auto world_pos = matrix_vector_multiply(mo->cs, pos);
    return jtk::vec3<float>(world_pos[0], world_pos[1], world_pos[2]);
    }
  return invalid_vertex;
  }

uint32_t view::get_index(int x, int y)
  {
  std::scoped_lock lock(_mut);
  if (x < 0 || y < 0 || x >= (int)_pixels.width() || y >= (int)_pixels.height())
    return (uint32_t)(-1);
  const auto& p = _pixels(x, y);
  if (p.db_id == 0)
    return (uint32_t)(-1);
  uint32_t closest_v = get_closest_vertex(p, get_vertices(_db, p.db_id), get_triangles(_db, p.db_id));    
  return closest_v;
  }

uint32_t view::get_id(int x, int y)
  {
  std::scoped_lock lock(_mut);
  if (x < 0 || y < 0 || x >= (int)_pixels.width() || y >= (int)_pixels.height())
    return (uint32_t)0;
  const auto& p = _pixels(x, y);
  if (p.db_id == 0)
    return (uint32_t)0;
  return p.db_id;
  }

void view::unzoom()
  {
  std::scoped_lock lock(_mut);
  ::unzoom(_scene);
  _refresh = true;
  }

uint32_t view::_get_semirandom_matcap_id(uint32_t object_id) const
  {
  return _matcap.get_semirandom_matcap_id_for_given_db_id(object_id);
  }

int64_t view::parametric(const std::array<double, 6>& domain, jtk::vec3<double>(*fun_ptr)(double, double))
  {
  std::scoped_lock lock(_mut);
  mesh* db_mesh;
  uint32_t id;
  _db.create_mesh(db_mesh, id);

  double vrange = domain[4] - domain[3];
  double urange = domain[1] - domain[0];

  int vsteps = std::round(vrange / domain[5]) + 1;
  int usteps = std::round(urange / domain[2]) + 1;

  db_mesh->vertices.resize((uint32_t)usteps*(uint32_t)vsteps);

  jtk::parallel_for((int)0, vsteps, [&](int v)
    {
    double vpar = domain[3] + domain[5] * v;
    for (int u = 0; u < usteps; ++u)
      {
      double upar = domain[0] + domain[2] * u;

      uint32_t idx = (uint32_t)v * (uint32_t)usteps + (uint32_t)u;
      auto res = fun_ptr(upar, vpar);
      db_mesh->vertices[idx][0] = (float)res[0];
      db_mesh->vertices[idx][1] = (float)res[1];
      db_mesh->vertices[idx][2] = (float)res[2];
      }
    });

  db_mesh->triangles.reserve((uint32_t)(usteps - 1)*(uint32_t)(vsteps - 1) * 2);
  db_mesh->uv_coordinates.reserve((uint32_t)(usteps - 1)*(uint32_t)(vsteps - 1) * 2);
  
  for (int v = 0; v < vsteps-1; ++v)  
    {
    double vpar = domain[3] + domain[5] * v;
    for (int u = 0; u < usteps-1; ++u)
      {
      double upar = domain[0] + domain[2] * u;
      uint32_t idx_00 = (uint32_t)v * (uint32_t)usteps + (uint32_t)u;
      uint32_t idx_01 = (uint32_t)(v+1) * (uint32_t)usteps + (uint32_t)u;
      uint32_t idx_11 = (uint32_t)(v + 1) * (uint32_t)usteps + (uint32_t)(u+1);
      uint32_t idx_10 = (uint32_t)v * (uint32_t)usteps + (uint32_t)(u + 1);
      db_mesh->triangles.emplace_back(idx_00, idx_10, idx_11);
      db_mesh->triangles.emplace_back(idx_00, idx_11, idx_01);

      jtk::vec3<jtk::vec2<float>> uv;
      uv[0][0] = (float)upar;
      uv[0][1] = (float)vpar;
      uv[1][0] = (float)(upar + domain[2]);
      uv[1][1] = (float)vpar;
      uv[2][0] = (float)(upar + domain[2]);
      uv[2][1] = (float)(vpar + domain[5]);
      db_mesh->uv_coordinates.push_back(uv);
      uv[1][1] = (float)(vpar + domain[5]);
      uv[2][0] = (float)upar;
      db_mesh->uv_coordinates.push_back(uv);
      }
    }

  for (auto& uv: db_mesh->uv_coordinates)
    {
    for (int i = 0; i < 3; ++i)
      {
      uv[i][0] -= domain[0];
      uv[i][1] -= domain[3];
      uv[i][0] /= urange;
      uv[i][1] /= vrange;
      }
    }

  db_mesh->texture = make_dummy_texture(usteps, vsteps, 1);

  db_mesh->cs = get_identity();
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  db_mesh->visible = true;
  add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return (int64_t)id;
  }

int64_t view::marching_cubes(const jtk::boundingbox3d<float>& bb, uint64_t width, uint64_t height, uint64_t depth, float isovalue, double(*fun_ptr)(double, double, double))
  {
  std::scoped_lock lock(_mut);
  mesh* db_mesh;
  uint32_t id;
  _db.create_mesh(db_mesh, id);
  jtk::marching_cubes(db_mesh->vertices, db_mesh->triangles, bb, (uint32_t)width, (uint32_t)height, (uint32_t)depth, isovalue, [&](double x, double y, double z)
    {
    return fun_ptr(x, y, z);
    }, [](float) {return true; });
  db_mesh->cs = get_identity();
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  db_mesh->visible = true;
  add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return (int64_t)id;
  }

void view::pointcloud_estimate_normals(uint32_t id, uint32_t k)
  {
  std::scoped_lock lock(_mut);
  pc* p = _db.get_pc(id);
  if (p)
    {
    p->normals = estimate_normals(*p, k);
    remove_object(id, _scene);
    if (p->visible)
      add_object(id, _scene, _db);
    _refresh = true;
    }
  }
  
void view::get_image(std::vector<uint32_t>& out, uint32_t& w, uint32_t& h, uint32_t id)
  {
  std::scoped_lock lock(_mut);
  out.clear();
  w = 0;
  h = 0;
  im* i = _db.get_image(id);
  if (i)
    {
    w = i->texture.width();
    h = i->texture.height();
    out.reserve(w*h);
    for (uint32_t y = 0; y < h; ++y)
      {
      uint32_t* p_im = i->texture.row(y);
      for (uint32_t x = 0; x < w; ++x, ++p_im)
        out.push_back(*p_im);
      }
    }
  }
  
int64_t view::image_pyr_down(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  im* i = _db.get_image(id);
  if (i)
    {
    im* new_object;
    uint32_t id_out;
    _db.create_image(new_object, id_out);
    new_object->texture = pyramid_down(i->texture);
    return id_out;
    }
  }
    
int64_t view::image_pyr_up(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  im* i = _db.get_image(id);
  if (i)
    {
    im* new_object;
    uint32_t id_out;
    _db.create_image(new_object, id_out);
    new_object->texture = pyramid_up(i->texture);
    return id_out;
    }
  return -1;
  }

int64_t view::image_gauss(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  im* i = _db.get_image(id);
  if (i)
    {
    im* new_object;
    uint32_t id_out;
    _db.create_image(new_object, id_out);
    new_object->texture = gauss(i->texture);
    return id_out;
    }
  return -1;
  }

bool view::write(uint32_t id, const char* filename)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (m)
    {
    return write_to_file(*m, filename);
    }
  pc* p = _db.get_pc((uint32_t)id);
  if (p)
    {
    return write_to_file(*p, filename);
    }
  mm* mo = _db.get_mm((uint32_t)id);
  if (mo)
    {
    return write_to_file(*mo, filename);
    }
  im* image = _db.get_image((uint32_t)id);
  if (image)
    {
    return write_to_file(*image, filename);
    }
  return false;
  }

std::vector<jtk::vec3<uint32_t>> view::triangles(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (m)
    return m->triangles;
  mm* mo = _db.get_mm((uint32_t)id);
  if (mo)
    return mo->shape.triangles;
  return std::vector<jtk::vec3<uint32_t>>();
  }

int64_t view::lscm(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh(id);
  if (m)
    {
    mesh* new_object;
    uint32_t new_id;
    _db.create_mesh(new_object, new_id);
    _matcap.map_db_id_to_matcap_id(new_id, _get_semirandom_matcap_id(new_id));
    *new_object = *m;
    new_object->visible = true;    
    
    auto lscm_uv = ::lscm(new_object->triangles.data(), (uint32_t)new_object->triangles.size(), new_object->vertices.data(), (uint32_t)new_object->vertices.size());
    optimize_aabb(lscm_uv, new_object->triangles.data(), (uint32_t)new_object->triangles.size(), (uint32_t)new_object->vertices.size());
    scale_to_unit(lscm_uv);
    if (new_object->texture.width() == 0 || new_object->texture.height() == 0)
      new_object->texture = make_dummy_texture(512, 512);

    new_object->uv_coordinates.clear();
    new_object->uv_coordinates.reserve(new_object->triangles.size());
    for (uint32_t t = 0; t < (uint32_t)new_object->triangles.size(); ++t)
      {
      jtk::vec3<jtk::vec2<float>> uv;
      for (uint32_t j = 0; j < 3; ++j)
        {
        const auto& vertex_uv = lscm_uv[new_object->triangles[t][j]];
        uv[j][0] = vertex_uv.first;
        uv[j][1] = vertex_uv.second;
        }
      new_object->uv_coordinates.push_back(uv);
      }

    if (new_object->visible)
      add_object(new_id, _scene, _db);
    prepare_scene(_scene);
    _refresh = true;
    return new_id;
    }
  return -1;
  }

int64_t view::make_minimal(const std::vector<jtk::vec3<float>>& vertices, uint32_t number_of_rings, uint32_t number_of_iterations)
  {
  std::scoped_lock lock(_mut);
  if (vertices.size() < 3)
    return -1;
  mesh* new_object;
  uint32_t new_id;
  _db.create_mesh(new_object, new_id);
  _matcap.map_db_id_to_matcap_id(new_id, _get_semirandom_matcap_id(new_id));
  new_object->vertices = vertices;
  new_object->cs = jtk::get_identity();
  new_object->visible = true;
  std::vector<uint32_t> hole;
  hole.reserve(vertices.size());
  for (uint32_t i = 0; i < (uint32_t)vertices.size(); ++i)
    hole.push_back(i);
  fill_hole_minimal_surface_parameters pars;
  pars.number_of_rings = number_of_rings;
  pars.iterations = number_of_iterations;
  fill_hole_minimal_surface(new_object->triangles, new_object->vertices, hole, pars);
  if (new_object->visible)
    add_object(new_id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return new_id;
  }

int64_t view::fill_hole(uint32_t id, const std::vector<uint32_t>& hole)
  {
  std::scoped_lock lock(_mut);
  std::vector<vec3<uint32_t>>* p_triangles = get_triangles(_db, id);
  std::vector<vec3<float>>* p_vertices = get_vertices(_db, id);
  if (p_triangles && p_vertices)
    {
    mesh* new_object;
    uint32_t new_id;
    _db.create_mesh(new_object, new_id);
    _matcap.map_db_id_to_matcap_id(new_id, _get_semirandom_matcap_id(new_id));
    new_object->triangles = *p_triangles;
    new_object->vertices = *p_vertices;
    new_object->visible = true;
    if (get_cs(_db, id))
      new_object->cs = *get_cs(_db, id);
    else
      new_object->cs = jtk::get_identity();

    jtk::mutable_adjacency_list adj_list((uint32_t)new_object->vertices.size(), new_object->triangles.data(), (uint32_t)new_object->triangles.size());
    jtk::fill_hole_ear<jtk::minimum_weight_ear>(new_object->triangles, new_object->vertices.data(), adj_list, hole);
    if (new_object->visible)
      add_object(new_id, _scene, _db);
    prepare_scene(_scene);
    _refresh = true;
    return new_id;
    }
  return -1;
  }

int64_t view::fill_hole_minimal(uint32_t id, const std::vector<uint32_t>& hole, uint32_t number_of_rings, uint32_t number_of_iterations)
  {
  std::scoped_lock lock(_mut);
  std::vector<vec3<uint32_t>>* p_triangles = get_triangles(_db, id);
  std::vector<vec3<float>>* p_vertices = get_vertices(_db, id);
  if (p_triangles && p_vertices)
    {
    mesh* new_object;
    uint32_t new_id;
    _db.create_mesh(new_object, new_id);
    _matcap.map_db_id_to_matcap_id(new_id, _get_semirandom_matcap_id(new_id));
    new_object->triangles = *p_triangles;
    new_object->vertices = *p_vertices;
    new_object->visible = true;
    if (get_cs(_db, id))
      new_object->cs = *get_cs(_db, id);
    else
      new_object->cs = jtk::get_identity();
    fill_hole_minimal_surface_parameters pars;
    pars.number_of_rings = number_of_rings;
    pars.iterations = number_of_iterations;
    fill_hole_minimal_surface(new_object->triangles, new_object->vertices, hole, pars);
    
    if (new_object->visible)
      add_object(new_id, _scene, _db);
    prepare_scene(_scene);
    _refresh = true;
    return new_id;
    }
  return -1;
  }

std::vector<std::vector<uint32_t>> view::holes(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  std::vector<std::vector<uint32_t>> holes_found;
  std::vector<vec3<uint32_t>>* p_triangles = get_triangles(_db, id);
  std::vector<vec3<float>>* p_vertices = get_vertices(_db, id);
  if (p_triangles && p_vertices)
    {
    jtk::adjacency_list adj_list((uint32_t)p_vertices->size(), p_triangles->data(), (uint32_t)p_triangles->size());
    holes_found = jtk::holes(adj_list, p_triangles->data());
    }
  return holes_found;
  }

void view::diagnose(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  std::vector<jtk::vec3<float>>* v = get_vertices(_db, id);
  std::vector<jtk::vec3<uint32_t>>* t = get_triangles(_db, id);
  if (v && t)
    {
    cork_options ops;
    ops.p_str = nullptr;
    ops.resolve_all_intersections = true;
    ops.use_parallel = true;
    ops.debug_folder = nullptr;
    ops.center = true;
    auto d = ::diagnose(*t, *v, ops);
    std::cout << "Diagnostics" << std::endl;
    std::cout << "-----------" << std::endl;
    std::cout << "number of triangles: " << d.number_of_triangles << std::endl;
    std::cout << "number of vertices: " << d.number_of_vertices << std::endl;
    std::cout << "number of self intersections: " << d.intersections << std::endl;
    std::cout << "number of degenerates: " << d.degenerate << std::endl;
    std::cout << "-----------" << std::endl;
    }
  }

int64_t view::resolve_intersections(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  std::vector<jtk::vec3<float>>* v = get_vertices(_db, id);
  std::vector<jtk::vec3<uint32_t>>* t = get_triangles(_db, id);
  if (v && t)
    {
    mesh* db_mesh;
    uint32_t id_out;
    _db.create_mesh(db_mesh, id_out);
    db_mesh->cs = *get_cs(_db, id);
    cork_options ops;
    ops.p_str = nullptr;
    ops.resolve_all_intersections = true;
    ops.use_parallel = true;
    ops.debug_folder = nullptr;
    ops.center = true;
    ::resolve_intersections(db_mesh->triangles, db_mesh->vertices, *t, *v, ops);
    db_mesh->visible = true;
    _matcap.map_db_id_to_matcap_id(id_out, _get_semirandom_matcap_id(id_out));
    if (db_mesh->visible)
      add_object(id_out, _scene, _db);
    //prepare_scene(_scene);
    //::unzoom(_scene);
    _refresh = true;
    return (int64_t)id_out;
    }
  return -1;
  }

int64_t view::csg(const std::vector<uint32_t>& ids, int csg_type)
  {
  std::scoped_lock lock(_mut);
  if (ids.size() < 2)
    return -1;
  uint32_t id1 = ids[0];
  uint32_t id2 = ids[1];
  std::vector<jtk::vec3<float>>* v1 = get_vertices(_db, id1);
  std::vector<jtk::vec3<float>>* v2 = get_vertices(_db, id2);
  std::vector<jtk::vec3<uint32_t>>* t1 = get_triangles(_db, id1);
  std::vector<jtk::vec3<uint32_t>>* t2 = get_triangles(_db, id2);
  jtk::float4x4* cs1 = get_cs(_db, id1);
  jtk::float4x4* cs2 = get_cs(_db, id2);

  if (v1 && v2 && t1 && t2 && cs1 && cs2)
    {
    mesh* db_mesh;
    uint32_t id;
    _db.create_mesh(db_mesh, id);
    db_mesh->cs = jtk::get_identity();
    cork_options ops;
    ops.p_str = nullptr;
    ops.resolve_all_intersections = false;
    ops.use_parallel = true;
    ops.debug_folder = nullptr;
    ops.center = true;

    switch (csg_type)
      {
      case 0:
      {
      compute_union(db_mesh->triangles, db_mesh->vertices, *t1, *v1, &(*cs1)[0], *t2, *v2, &(*cs2)[0], ops);
      break;
      }
      case 1:
      {
      compute_difference(db_mesh->triangles, db_mesh->vertices, *t1, *v1, &(*cs1)[0], *t2, *v2, &(*cs2)[0], ops);
      break;
      }
      case 2:
      {
      compute_intersection(db_mesh->triangles, db_mesh->vertices, *t1, *v1, &(*cs1)[0], *t2, *v2, &(*cs2)[0], ops);
      break;
      }
      }

    for (uint32_t j = 2; j < ids.size(); ++j)
      {
      v2 = get_vertices(_db, ids[j]);
      t2 = get_triangles(_db, ids[j]);
      cs2 = get_cs(_db, ids[j]);
      if (v2 && t2 && cs2)
        {
        mesh tmp;
        tmp.vertices.swap(db_mesh->vertices);
        tmp.triangles.swap(db_mesh->triangles);

        switch (csg_type)
          {
          case 0:
          {
          compute_union(db_mesh->triangles, db_mesh->vertices, tmp.triangles, tmp.vertices, &(db_mesh->cs)[0], *t2, *v2, &(*cs2)[0], ops);
          break;
          }
          case 1:
          {
          compute_difference(db_mesh->triangles, db_mesh->vertices, tmp.triangles, tmp.vertices, &(db_mesh->cs)[0], *t2, *v2, &(*cs2)[0], ops);
          break;
          }
          case 2:
          {
          compute_intersection(db_mesh->triangles, db_mesh->vertices, tmp.triangles, tmp.vertices, &(db_mesh->cs)[0], *t2, *v2, &(*cs2)[0], ops);
          break;
          }
          }
        }
      }
    
    db_mesh->visible = true;
    _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
    if (db_mesh->visible)
      add_object(id, _scene, _db);
    //prepare_scene(_scene);
    //::unzoom(_scene);
    _refresh = true;
    return (int64_t)id;
    }
  return -1;
  }

std::vector<jtk::vec3<float>> view::vertexnormals(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  pc* p = _db.get_pc((uint32_t)id);
  if (p)
    return p->normals;
  mesh* m = _db.get_mesh((uint32_t)id);
  if (m)
    {
    std::vector<vec3<float>> triangle_normals, vertex_normals;
    compute_triangle_normals(triangle_normals, m->vertices.data(), m->triangles.data(), m->triangles.size());
    compute_vertex_normals(vertex_normals, triangle_normals.data(), m->vertices.data(), m->vertices.size(), m->triangles.data(), m->triangles.size());
    return vertex_normals;
    }
  mm* m2 = _db.get_mm((uint32_t)id);
  if (m2)
    {
    std::vector<vec3<float>> triangle_normals, vertex_normals;
    compute_triangle_normals(triangle_normals, m2->vertices.data(), m2->shape.triangles.data(), m2->shape.triangles.size());
    compute_vertex_normals(vertex_normals, triangle_normals.data(), m2->vertices.data(), m2->vertices.size(), m2->shape.triangles.data(), m2->shape.triangles.size());
    return vertex_normals;
    }
  return std::vector<jtk::vec3<float>>();
  }

std::vector<jtk::vec3<float>> view::trianglenormals(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (m)
    {
    std::vector<vec3<float>> triangle_normals;
    compute_triangle_normals(triangle_normals, m->vertices.data(), m->triangles.data(), m->triangles.size());    
    return triangle_normals;
    }
  mm* m2 = _db.get_mm((uint32_t)id);
  if (m2)
    {
    std::vector<vec3<float>> triangle_normals;
    compute_triangle_normals(triangle_normals, m2->vertices.data(), m2->shape.triangles.data(), m2->shape.triangles.size());    
    return triangle_normals;
    }
  return std::vector<jtk::vec3<float>>();
  }

std::vector<uint32_t> view::vertexcolors(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (m)
    return convert_vertex_colors(m->vertex_colors);
  pc* p = _db.get_pc((uint32_t)id);
  if (p)
    return p->vertex_colors;
  return std::vector<uint32_t>();
  }

std::vector<jtk::vec3<float>> view::vertices(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (m)
    return m->vertices;
  pc* p = _db.get_pc((uint32_t)id);
  if (p)
    return p->vertices;
  mm* mo = _db.get_mm((uint32_t)id);
  if (mo)
    return mo->vertices;
  return std::vector<jtk::vec3<float>>();
  }

int64_t view::mm_coeff_size(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return -1;
  return (int64_t)m->shape.U.cols();
  }

int64_t view::mm_shape_size(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return -1;
  return (int64_t)m->shape.U.rows();
  }

double view::mm_sigma(uint32_t id, int64_t idx)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return std::numeric_limits<double>::quiet_NaN();
  return jtk::sigma(m->shape, idx);
  }

std::vector<float> view::mm_coeff(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return std::vector<float>();
  return m->coefficients;
  }

std::vector<float> view::mm_basic_shape_coeff(uint32_t id, int64_t shape_id)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return std::vector<float>();
  return jtk::get_basic_shape(m->shape, (uint32_t)shape_id);
  }

void view::mm_coeff_set(uint32_t id, const std::vector<float>& coeff)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return;
  m->coefficients = coeff;
  m->vertices = jtk::get_vertices(m->shape, m->coefficients);
  remove_object(id, _scene);
  if (m->visible)
    {
    add_object(id, _scene, _db);
    _refresh = true;
    }
  }

int64_t view::mm_color_coeff_size(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return -1;
  return (int64_t)m->color.U.cols();
  }

int64_t view::mm_color_shape_size(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return -1;
  return (int64_t)m->color.U.rows();
  }

double view::mm_color_sigma(uint32_t id, int64_t idx)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return std::numeric_limits<double>::quiet_NaN();
  return jtk::sigma(m->color, idx);
  }

std::vector<float> view::mm_color_coeff(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return std::vector<float>();
  return m->color_coefficients;
  }

std::vector<float> view::mm_color_basic_shape_coeff(uint32_t id, int64_t shape_id)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return std::vector<float>();
  return jtk::get_basic_shape(m->color, (uint32_t)shape_id);
  }

void view::mm_color_coeff_set(uint32_t id, const std::vector<float>& coeff)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)id);
  if (!m)
    return;
  m->color_coefficients = coeff;
  m->vertex_colors = jtk::get_vertices(m->color, m->color_coefficients);
  clamp_vertex_colors(m->vertex_colors);
  remove_object(id, _scene);
  if (m->visible)
    {
    add_object(id, _scene, _db);
    _refresh = true;
    }
  }

int64_t view::mm_to_mesh(int32_t mm_id)
  {
  std::scoped_lock lock(_mut);
  mm* m = _db.get_mm((uint32_t)mm_id);
  if (!m)
    return -1;
  mesh* db_mesh;
  uint32_t id;
  _db.create_mesh(db_mesh, id);
  db_mesh->vertices = m->vertices;
  db_mesh->triangles = m->shape.triangles;
  db_mesh->vertex_colors = m->vertex_colors;
  db_mesh->cs = m->cs;
  db_mesh->visible = true;
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  if (db_mesh->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return (int64_t)id;
  }

std::vector<jtk::vec3<uint8_t>> view::mesh_texture_to_vertexcolors(uint32_t id)
  {
  std::vector<jtk::vec3<uint8_t>> colors;
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (m && !m->uv_coordinates.empty())
    {
    colors.resize(m->vertices.size());
    for (uint32_t t = 0; t < (uint32_t)m->triangles.size(); ++t)
      {
      const auto& tria = m->triangles[t];
      const auto& uv = m->uv_coordinates[t];
      for (int j = 0; j < 3; ++j)
        {
        float u = uv[j][0];
        float v = uv[j][1];
        if (u >= 0.f && u <= 1.f && v >= 0.f && v <= 1.f)
          {
          int U = (int)((m->texture.width()-1)*u);
          int V = (int)((m->texture.height() - 1)*v);
          uint32_t clr = m->texture(U, V);
          colors[tria[j]][0] = clr & 255;
          colors[tria[j]][1] = (clr >> 8) & 255;
          colors[tria[j]][2] = (clr >> 16) & 255;
          }
        }
      }
    }
  return colors;
  }

void view::shape_predictor_set_flip_horizontal(uint32_t id, bool flip)
  {
  std::scoped_lock lock(_mut);
  sp* m = _db.get_sp((uint32_t)id);
  if (m)
    {
    m->flip_horizontal = flip;
    _refresh = true;
    }
  }

void view::shape_predictor_link(uint32_t id, sp::object_detector_link link)
  {
  std::scoped_lock lock(_mut);
  sp* m = _db.get_sp((uint32_t)id);
  if (m)
    {
    m->odl = link;
    }
  }

int64_t view::load_shape_predictor(const char* filename)
  {
  std::scoped_lock lock(_mut);
  sp m;
  std::string f(filename);
  bool res = read_from_file(m, f);
  if (!res)
    return -1;
  sp* db_sp;
  uint32_t id;
  _db.create_sp(db_sp, id);
  swap(*db_sp, m);
  return (int64_t)id;
  }

void view::set_show_ear_left_detector(bool b)
  {
  std::scoped_lock lock(_mut);
  _show_ear_left_detector = b;
  _refresh = true;
  }

void view::set_show_ear_right_detector(bool b)
  {
  std::scoped_lock lock(_mut);
  _show_ear_right_detector = b;
  _refresh = true;
  }

void view::set_show_face_detector(bool b)
  {
  std::scoped_lock lock(_mut);
  _show_face_detector = b;
  _refresh = true;
  }

std::vector<std::pair<long, long>> view::shape_predict(uint32_t id, const rect& r)
  {
  std::scoped_lock lock(_mut);
  std::vector<std::pair<long, long>> landmarks;
  bool show_left_ear = _show_ear_left_detector;
  bool show_right_ear = _show_ear_right_detector;
  bool show_face = _show_face_detector;
  sp* m = _db.get_sp((uint32_t)id);
  if (!m)
    return landmarks;
  if (m->p_shape_predictor.get())
    {
    _show_ear_left_detector = false;
    _show_ear_right_detector = false;
    _show_face_detector = false;
    render_scene();
    landmarks = predict(*m, r, _canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data());
    }
  _show_ear_left_detector = show_left_ear;
  _show_ear_right_detector = show_right_ear;
  _show_face_detector = show_face;
  _refresh = _show_ear_left_detector || _show_ear_right_detector || _show_face_detector;
  return landmarks;
  }

std::vector<rect> view::ear_detect(bool right_ear)
  {
  std::scoped_lock lock(_mut);
  std::vector<rect> rectangles;
  bool show_left_ear = _show_ear_left_detector;
  bool show_right_ear = _show_ear_right_detector;
  bool show_face = _show_face_detector;

  if (p_ear_detector.get())
    {
    _show_ear_left_detector = false;
    _show_ear_right_detector = false;
    _show_face_detector = false;
    render_scene();
    rectangles = p_ear_detector->detect(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data(), right_ear);
    }
  _show_ear_left_detector = show_left_ear;
  _show_ear_right_detector = show_right_ear;
  _show_face_detector = show_face;
  _refresh = _show_ear_left_detector || _show_ear_right_detector || _show_face_detector;
  return rectangles;
  }

std::vector<rect> view::face_detect()
  {
  std::scoped_lock lock(_mut);
  std::vector<rect> rectangles;
  bool show_left_ear = _show_ear_left_detector;
  bool show_right_ear = _show_ear_right_detector;
  bool show_face = _show_face_detector;
  if (p_face_detector.get())
    {
    _show_ear_left_detector = false;
    _show_ear_right_detector = false;
    _show_face_detector = false;
    render_scene();
    rectangles = p_face_detector->detect(_canvas.get_image().width(), _canvas.get_image().height(), _canvas.get_image().stride(), _canvas.get_image().data());
    }
  _show_ear_left_detector = show_left_ear;
  _show_ear_right_detector = show_right_ear;
  _show_face_detector = show_face;
  _refresh = _show_ear_left_detector || _show_ear_right_detector || _show_face_detector;
  return rectangles;
  }

void view::fit_mm_to_partial_positions(uint32_t mm_id, const std::vector<uint32_t>& vertex_indices, const std::vector<jtk::vec3<float>>& vertex_positions)
  {
  std::scoped_lock lock(_mut);
  mm* morph = _db.get_mm((uint32_t)mm_id);
  if (!morph)
    return;
  fit_to_partial_positions(*morph, vertex_indices, vertex_positions);
  remove_object(mm_id, _scene);
  if (morph->visible)
    {
    add_object(mm_id, _scene, _db);
    _refresh = true;
    }
  }

void view::fit_mm(uint32_t mm_id, uint32_t mesh_id)
  {
  std::scoped_lock lock(_mut); 
  mesh* m = _db.get_mesh((uint32_t)mesh_id);
  if (!m)
    return;
  mm* morph = _db.get_mm((uint32_t)mm_id);
  if (!morph)
    return;
  fit_to_mesh(*morph, *m);
  remove_object(mm_id, _scene);
  if (morph->visible)
    {
    add_object(mm_id, _scene, _db);
    _refresh = true;
    }
  }

int64_t view::poisson(uint32_t pc_id, uint32_t depth)
  {
  std::scoped_lock lock(_mut);
  pc* p = _db.get_pc(pc_id);
  if (!p)
    return -1;  
  if (p->normals.empty())
    {
    std::cout << "error: This pointcloud has no normals\n";
    return -1;
    }
  poisson_reconstruction_screened_parameters pars;
  pars.depth = depth;  
  mesh* db_mesh;
  uint32_t id;
  _db.create_mesh(db_mesh, id);
  if (p->vertex_colors.empty())  
    poisson_reconstruction_screened(db_mesh->vertices, db_mesh->triangles, p->vertices, p->normals, pars);
  else
    poisson_reconstruction_screened(db_mesh->vertices, db_mesh->triangles, db_mesh->vertex_colors, p->vertices, p->normals, p->vertex_colors, pars);
  db_mesh->cs = jtk::get_identity();
  db_mesh->visible = true;
  _matcap.map_db_id_to_matcap_id(id, _get_semirandom_matcap_id(id));
  if (db_mesh->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
  return (int64_t)id;
  }

void view::info(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  pc* p = _db.get_pc(id);
  if (p)
    ::info(*p);
  mesh* m = _db.get_mesh(id);
  if (m)
    ::info(*m);
  mm* morph = _db.get_mm(id);
  if (morph)
    ::info(*morph);
  sp* shape_pred = _db.get_sp(id);
  if (shape_pred)
    ::info(*shape_pred);
  im* image = _db.get_image(id);
  if (image)
    ::info(*image);
  }

void view::cs_apply(uint32_t id)
  {
  std::scoped_lock lock(_mut);
  pc* p = _db.get_pc(id);
  bool visible = false;
  if (p)
    {
    ::cs_apply(*p);
    visible = p->visible;
    }
  mesh* m = _db.get_mesh(id);
  if (m)
    {
    ::cs_apply(*m);
    visible = m->visible;
    }
  mm* morph = _db.get_mm(id);
  if (morph)
    {
    ::cs_apply(*morph);
    visible = morph->visible;
    }

  remove_object(id, _scene);
  if (visible)
    {
    add_object(id, _scene, _db);
    _refresh = true;
    }
  }

void view::poll_for_events()
  {
  _m.right_button_down = false;
  _m.left_button_down = false;
  _m.wheel_down = false;
  SDL_Event event;
  while (SDL_PollEvent(&event))
    {
    if (event.type == SDL_QUIT)
      {
      SDL_FlushEvents(SDL_FIRSTEVENT, SDL_LASTEVENT);
      this->hide();
      }
    if (event.type == SDL_MOUSEMOTION)
      {
      _m.prev_mouse_x = _m.mouse_x;
      _m.prev_mouse_y = _m.mouse_y;
      _m.mouse_x = float(event.motion.x);
      _m.mouse_y = float(event.motion.y);
      }
    if (event.type == SDL_MOUSEBUTTONDOWN)
      {
      if (event.button.button == 2)
        {
        _m.wheel_mouse_pressed = true;
        _m.wheel_down = true;
        }
      else if (event.button.button == 1)
        {
        _m.left_dragging = true;
        _m.left_button_down = true;
        }
      else if (event.button.button == 3)
        {
        _m.right_dragging = true;
        _m.right_button_down = true;
        }
      }
    if (event.type == SDL_MOUSEBUTTONUP)
      {
      if (event.button.button == 2)
        _m.wheel_mouse_pressed = false;
      else if (event.button.button == 1)
        {
        _m.left_dragging = false;
        }
      else if (event.button.button == 3)
        {
        _m.right_dragging = false;
        }
      }
    if (event.type == SDL_MOUSEWHEEL)
      {
      _m.wheel_rotation += event.wheel.y;
      }
    if (event.type == SDL_KEYDOWN)
      {
      switch (event.key.keysym.sym)
        {
        case SDLK_LCTRL:
        case SDLK_RCTRL:
        {
        _m.ctrl_pressed = true;
        break;
        }
        }
      }
    if (event.type == SDL_KEYUP)
      {
      switch (event.key.keysym.sym)
        {
        case SDLK_ESCAPE:
        {
        _quit = true;
        break;
        }
        case SDLK_LCTRL:
        case SDLK_RCTRL:
        {
        _m.ctrl_pressed = false;
        break;
        }
        case SDLK_1:
        {
        resize_canvas(512, 512);
        break;
        }
        case SDLK_2:
        {
        resize_canvas(800, 600);
        break;
        }
        case SDLK_3:
        {
        resize_canvas(1024, 768);
        break;
        }
        case SDLK_4:
        {
        resize_canvas(1600, 900);
        break;
        }
        case SDLK_s:
        {
        _settings.shadow = !_settings.shadow;
        _refresh = true;
        break;
        }
        case SDLK_t:
        {
        _settings.textured = !_settings.textured;
        _refresh = true;
        break;
        }
        case SDLK_l:
        {
        _settings.shading = !_settings.shading;
        _refresh = true;
        break;
        }
        case SDLK_e:
        {
        _settings.edges = !_settings.edges;
        _refresh = true;
        break;
        }
        case SDLK_b:
        {
        _settings.one_bit = !_settings.one_bit;
        _refresh = true;
        break;
        }
        case SDLK_u:
        {
        ::unzoom(_scene);
        _refresh = true;
        break;
        }
        case SDLK_v:
        {
        _settings.vertexcolors = !_settings.vertexcolors;
        _refresh = true;
        break;
        }
        case SDLK_w:
        {
        _settings.wireframe = !_settings.wireframe;
        _refresh = true;
        break;
        }
        case SDLK_SPACE:
        {
        if (mouse_in_canvas())
          {
          pixel p_actual;
          _canvas.get_pixel(p_actual, _m, (float)_canvas_pos_x, (float)_canvas_pos_y);
          //if (p_actual.object_id != (uint32_t)(-1))
          if (p_actual.db_id)
            {
            uint32_t closest_v = get_closest_vertex(p_actual, get_vertices(_db, p_actual.db_id), get_triangles(_db, p_actual.db_id));
            auto V = (*get_vertices(_db, p_actual.db_id))[closest_v];
            std::cout << std::endl;
            std::cout << "---------------------------------------" << std::endl;
            std::cout << "database id: " << p_actual.db_id << std::endl;
            std::cout << "vertex index: " << closest_v << std::endl;
            std::cout << "vertex coordinates: " << V[0] << ", " << V[1] << ", " << V[2] << ")" << std::endl;
            std::cout << "---------------------------------------" << std::endl;
            }
          }
        break;
        }
        }
      }
    if (event.type == SDL_WINDOWEVENT)
      {
      if (event.window.event == SDL_WINDOWEVENT_SIZE_CHANGED)
        {
        int new_w, new_h;
        SDL_GetWindowSize(_window, &new_w, &new_h);
        _w = (uint32_t)new_w < _w_max ? (uint32_t)new_w : _w_max;
        _h = (uint32_t)new_h < _h_max ? (uint32_t)new_h : _h_max;
        _screen = image<uint32_t>(_w, _h);
        clear_screen(_screen);
        _canvas_pos_x = ((int32_t)_w - (int32_t)_canvas.width()) / 2;
        if (_canvas_pos_x & 3)
          _canvas_pos_x += 4 - (_canvas_pos_x & 3);
        _canvas_pos_y = ((int32_t)_h - (int32_t)_canvas.height()) / 2;
        _refresh = true;
        }
      }
    }
  }

void view::blit_screen_to_opengl_texture()
  {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, _w, 0, _h, 0, 1);
  glMatrixMode(GL_MODELVIEW);
  glViewport(0, 0, _w, _h);
  glClear(GL_COLOR_BUFFER_BIT);

  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, _gl_texture);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);


  const uint32_t step_block_x = _w / _gl_texture_w + 1;
  const uint32_t step_block_y = _h / _gl_texture_h + 1;

  for (uint32_t block_y = 0; block_y < step_block_y; ++block_y)
    {
    for (uint32_t block_x = 0; block_x < step_block_x; ++block_x)
      {
      const uint32_t x_offset = block_x * _gl_texture_w;
      const uint32_t y_offset = block_y * _gl_texture_h;
      uint32_t wi = std::min<uint32_t>(_gl_texture_w, _w - x_offset);
      uint32_t he = std::min<uint32_t>(_gl_texture_h, _h - y_offset);

      uint32_t* p_im_start = (uint32_t*)(_screen.row(y_offset) + x_offset);
      glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
      glPixelStorei(GL_UNPACK_ROW_LENGTH, _screen.stride());
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, wi, he, GL_RGBA, GL_UNSIGNED_BYTE, (void*)p_im_start);
      gl_check_error("glTexSubImage2D in view.cpp");
      const float tex_w = float(wi) / float(_gl_texture_w);
      const float tex_h = float(he) / float(_gl_texture_h);

      const float x0 = float(x_offset);
      const float y0 = float(y_offset);

      const float x1 = x0 + float(wi);
      const float y1 = y0 + float(he);
      glBegin(GL_TRIANGLE_STRIP);
      glTexCoord2f(0.0, 0.0);
      glVertex2f(x0, y0);
      glTexCoord2f(tex_w, 0.0);
      glVertex2f(x1, y0);
      glTexCoord2f(0.0, tex_h);
      glVertex2f(x0, y1);
      glTexCoord2f(tex_w, tex_h);
      glVertex2f(x1, y1);
      glEnd();
      gl_check_error("GL_TRIANGLE_STRIP in view.cpp");
      }
    }
  }

void view::resize_canvas(uint32_t canvas_w, uint32_t canvas_h)
  {
  _refresh = true;
  if (canvas_w > _w)
    canvas_w = _w;
  if (canvas_h > _h)
    canvas_h = _h;
  _canvas.resize(canvas_w, canvas_h);
  _canvas_pos_x = ((int32_t)_w - (int32_t)_canvas.width()) / 2;
  if (_canvas_pos_x & 3)
    _canvas_pos_x += 4 - (_canvas_pos_x & 3);
  _canvas_pos_y = ((int32_t)_h - (int32_t)_canvas.height()) / 2;
  }

bool view::mouse_in_canvas() const
  {
  return ((_m.mouse_x >= (float)_canvas_pos_x) && (_m.mouse_y >= (float)_canvas_pos_y) &&
    (_m.mouse_x < (float)(_canvas_pos_x + _canvas.width())) && (_m.mouse_y < (float)(_canvas_pos_y + _canvas.height())));
  }

void view::do_canvas_mouse()
  {
  if (!mouse_in_canvas())
    return;
  _canvas.do_mouse(_refresh, _m, _scene, (float)_canvas_pos_x, (float)_canvas_pos_y);
  }

namespace
  {
  void _draw_square(int x, int y, int size, jtk::image<uint32_t>& im, uint32_t clr)
    {
    int hs = (size - 1) / 2;
    if (x >= hs && y >= hs && x < ((int)im.width() - hs) && y < ((int)im.height() - hs))
      {
      for (int a = x - hs; a <= x + hs; ++a)
        {
        im(a, y - hs) = clr;
        im(a, y + hs) = clr;
        }
      for (int a = y - hs; a <= y + hs; ++a)
        {
        im(x - hs, a) = clr;
        im(x + hs, a) = clr;
        }
      }
    }
  }

void view::render_mouse()
  {
  if (!mouse_in_canvas())
    return;
  pixel p_actual;
  const float pos_x = _m.mouse_x;
  const float pos_y = _m.mouse_y;
  _canvas.get_pixel(p_actual, _m, (float)_canvas_pos_x, (float)_canvas_pos_y);
  const uint32_t clr = 0xff0000d0;
  //if (p_actual.object_id != (uint32_t)(-1))
  if (p_actual.db_id)
    {
    uint32_t closest_v = get_closest_vertex(p_actual, get_vertices(_db, p_actual.db_id), get_triangles(_db, p_actual.db_id));
    auto V = (*get_vertices(_db, p_actual.db_id))[closest_v];
    jtk::float4 V4(V[0], V[1], V[2], 1.f);
    V4 = jtk::matrix_vector_multiply(*get_cs(_db, p_actual.db_id), V4);
    V4 = jtk::matrix_vector_multiply(_scene.coordinate_system_inv, V4);
    V4 = jtk::matrix_vector_multiply(_canvas.get_projection_matrix(), V4);
    int x = (int)(((V4[0] / V4[3] + 1.f) / 2.f)*_canvas.width() - 0.5f);
    int y = (int)(((V4[1] / V4[3] + 1.f) / 2.f)*_canvas.height() - 0.5f);
    _draw_square(x + _canvas_pos_x, _screen.height() - 1 - (y + _canvas_pos_y), 3, _screen, clr);
    _draw_square(x + _canvas_pos_x, _screen.height() - 1 - (y + _canvas_pos_y), 5, _screen, 0xff000000);
    }
  }

void view::quit()
  {
  _quit = true;
  }

void view::loop()
  {
  while (!_quit)
    {

    if (_window)
      {
      poll_for_events();
      do_canvas_mouse();
      }


    if (_refresh)
      {
          {
          std::scoped_lock lock(_mut);
          render_scene();
          }
      }

    clear_screen(_screen);
    _canvas.blit_onto(_screen, _canvas_pos_x, _canvas_pos_y);

    if (_window)
      {
      render_mouse();
      blit_screen_to_opengl_texture();
      SDL_GL_SwapWindow(_window);
      }
    else
      std::this_thread::sleep_for(std::chrono::duration<double, std::milli>(16.0));

    if (_suspend)
      {
      std::scoped_lock lock(_mut);
      _suspend = false;
      delete_window();
      }

    if (_resume)
      {
      std::scoped_lock lock(_mut);
      _resume = false;
      prepare_window();
      }
    }
  }
