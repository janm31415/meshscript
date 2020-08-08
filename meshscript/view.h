#pragma once

#include <SDL.h>
#include <SDL_opengl.h>

#include "canvas.h"
#include <jtk/image.h>
#include "matcap.h"
#include "db.h"
#include "scene.h"
#include "mouse.h"

#include <jtk/qbvh.h>

#include <mutex>

class view
  {
  public:
    view();
    ~view();

    int64_t load_mesh_from_file(const char* filename);

    int64_t load_mesh(const std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles);

    void set_coordinate_system(uint32_t id, const jtk::float4x4& cs);

    void premultiply_coordinate_system(uint32_t id, const jtk::float4x4& cs);

    jtk::float4x4 get_coordinate_system(uint32_t id);

    void set_vertex_colors(uint32_t id, const std::vector<jtk::vec3<uint8_t>>& colors);

    bool vertices_to_csv(int64_t id, const char* filename);

    bool triangles_to_csv(int64_t id, const char* filename);

    void hide(int64_t id);

    void show(int64_t id);

    void hide();

    void show();

    void set_color(int64_t id, int64_t clr_id);

    void set_shading(bool b);

    void set_edges(bool b);

    void set_textured(bool b);

    void set_shadow(bool b);

    void set_wireframe(bool b);

    void set_one_bit(bool b);

    void set_image_size(int w, int h);

    jtk::image<uint32_t> get_image();

    void set_bg_color(uint8_t r, uint8_t g, uint8_t b);

    jtk::vec3<float> get_world_position(int x, int y);

    void unzoom();

    void loop();

    void quit();

  private:

    void poll_for_events();

    void blit_screen_to_opengl_texture();

    void resize_canvas(uint32_t canvas_w, uint32_t canvas_h);

    bool mouse_in_canvas() const;

    void do_canvas_mouse();

    void prepare_window();

    void delete_window();

    void render_mouse();

  private:

    SDL_Window* _window;
    uint32_t _w, _h;
    uint32_t _w_max, _h_max;
    mouse_data _m;
    bool _quit;
    bool _refresh;
    bool _suspend;
    bool _resume;
    scene _scene;
    db _db;

    GLuint _gl_texture;
    uint32_t _gl_texture_w;
    uint32_t _gl_texture_h;

    jtk::image<uint32_t> _screen;
    canvas _canvas;
    canvas::canvas_settings _settings;
    int32_t _canvas_pos_x, _canvas_pos_y;

    jtk::image<pixel> _pixels;
    matcapmap _matcap;

    std::mutex _mut;
  };