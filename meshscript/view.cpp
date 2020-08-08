#include "view.h"
#include "mesh.h"
#include "pc.h"

#include <iostream>

#include <SDL_syswm.h>

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

  //prepare_window();

  _screen = image<uint32_t>(_w, _h);
  clear_screen(_screen);

  _m.left_dragging = false;
  _m.right_dragging = false;
  _m.right_button_down = false;
  _m.left_button_down = false;
  _m.wheel_down = false;
  _m.wheel_mouse_pressed = false;
  _m.mouse_x = 0.f;
  _m.mouse_y = 0.f;
  _m.prev_mouse_x = 0.f;
  _m.prev_mouse_y = 0.f;
  _m.wheel_rotation = 0.f;

  _canvas = canvas(800, 600);
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

  _suspend = false;
  _resume = false;
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
  _matcap.map_db_id_to_matcap[id] = (id % _matcap.matcaps.size());
  if (db_mesh->visible)
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
  _matcap.map_db_id_to_matcap[id] = (id % _matcap.matcaps.size());
  if (db_pc->visible)
    add_object(id, _scene, _db);
  prepare_scene(_scene);
  ::unzoom(_scene);
  _refresh = true;
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
  _matcap.map_db_id_to_matcap[id] = (id % _matcap.matcaps.size());
  db_mesh->visible = true;
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
  if (!m)
    return;
  m->cs = cs;
  remove_object(id, _scene);
  if (m->visible)
    add_object(id, _scene, _db);

  prepare_scene(_scene);
  _refresh = true;
  }

void view::render_scene()
  {
  // assumes a lock has been set already
  _canvas.update_settings(_settings);
  _canvas.render_scene(&_scene);
  _pixels = _canvas.get_pixels();
  _canvas.canvas_to_image(_pixels, _matcap);
  _refresh = false;
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
  if (!m)
    return jtk::get_identity();
  return m->cs;
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
  if (!m)
    return;
  m->cs = jtk::matrix_matrix_multiply(cs, m->cs);
  remove_object(id, _scene);
  if (m->visible)
    add_object(id, _scene, _db);

  prepare_scene(_scene);
  _refresh = true;
  }

void view::set_vertex_colors(uint32_t id, const std::vector<jtk::vec3<uint8_t>>& colors)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh(id);
  if (!m)
    return;
  m->vertex_colors.clear();
  m->vertex_colors.reserve(colors.size());
  for (auto& clr : colors)
    m->vertex_colors.emplace_back(clr[0] / 255.f, clr[1] / 255.f, clr[2] / 255.f);
  remove_object(id, _scene);
  if (m->visible)
    add_object(id, _scene, _db);
  _refresh = true;
  }

bool view::vertices_to_csv(int64_t id, const char* filename)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (!m)
    return false;
  std::string fn(filename);
  return ::vertices_to_csv(*m, fn);
  }

bool view::triangles_to_csv(int64_t id, const char* filename)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (!m)
    return false;
  std::string fn(filename);
  return ::triangles_to_csv(*m, fn);
  }

void view::hide(int64_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (!m)
    return;
  m->visible = false;
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

void view::set_color(int64_t id, int64_t clr_id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (!m)
    return;
  _matcap.map_db_id_to_matcap[(uint32_t)id] = (uint32_t)(clr_id % _matcap.matcaps.size());
  _refresh = true;
  }

void view::show(int64_t id)
  {
  std::scoped_lock lock(_mut);
  mesh* m = _db.get_mesh((uint32_t)id);
  if (!m)
    return;
  m->visible = true;
  remove_object((uint32_t)id, _scene);
  add_object((uint32_t)id, _scene, _db);
  prepare_scene(_scene);
  _refresh = true;
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

jtk::vec3<float> view::get_world_position(int x, int y)
  {
  std::scoped_lock lock(_mut);
  jtk::vec3<float> invalid_vertex(std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN());
  if (x < 0 || y < 0 || x >= (int)_pixels.width() || y >= (int)_pixels.height())
    return invalid_vertex;
  const auto& p = _pixels(x, y);
  if (p.triangle_id == (uint32_t)-1)
    return invalid_vertex;
  mesh* m = _db.get_mesh((uint32_t)p.db_id);
  if (!m)
    return invalid_vertex;
  const uint32_t v0 = m->triangles[p.triangle_id][0];
  const uint32_t v1 = m->triangles[p.triangle_id][1];
  const uint32_t v2 = m->triangles[p.triangle_id][2];
  const float4 V0(m->vertices[v0][0], m->vertices[v0][1], m->vertices[v0][2], 1.f);
  const float4 V1(m->vertices[v1][0], m->vertices[v1][1], m->vertices[v1][2], 1.f);
  const float4 V2(m->vertices[v2][0], m->vertices[v2][1], m->vertices[v2][2], 1.f);
  const float4 pos = V0 * (1.f - p.barycentric_u - p.barycentric_v) + p.barycentric_u*V1 + p.barycentric_v*V2;
  auto world_pos = matrix_vector_multiply(m->cs, pos);
  return jtk::vec3<float>(world_pos[0], world_pos[1], world_pos[2]);
  }

void view::unzoom()
  {
  std::scoped_lock lock(_mut);
  ::unzoom(_scene);
  _refresh = true;
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
      _quit = true;
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
    if (event.type == SDL_KEYUP)
      {
      switch (event.key.keysym.sym)
        {
        case SDLK_ESCAPE:
        {
        _quit = true;
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
          if (p_actual.triangle_id != (uint32_t)(-1))
            {
            uint32_t closest_v = get_closest_vertex(p_actual, get_vertices(_db, p_actual.db_id)->data(), get_triangles(_db, p_actual.db_id)->data());
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
  _canvas = canvas(canvas_w, canvas_h);
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
  if (p_actual.triangle_id != (uint32_t)(-1))
    {
    uint32_t closest_v = get_closest_vertex(p_actual, get_vertices(_db, p_actual.db_id)->data(), get_triangles(_db, p_actual.db_id)->data());
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