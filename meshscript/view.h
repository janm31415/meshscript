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
#include <memory>

class face_detector;

class view
  {
  public:
    view();
    ~view();

    int64_t load_mesh_from_file(const char* filename);

    int64_t load_morphable_model_from_file(const char* filename);

    int64_t load_pc_from_file(const char* filename);

    int64_t load_mesh(const std::vector<jtk::vec3<float>>& vertices, const std::vector<jtk::vec3<uint32_t>>& triangles);

    void set_coordinate_system(uint32_t id, const jtk::float4x4& cs);

    void set_coordinate_system(const jtk::float4x4& cs);

    void premultiply_coordinate_system(uint32_t id, const jtk::float4x4& cs);

    jtk::float4x4 get_coordinate_system(uint32_t id);

    jtk::float4x4 get_coordinate_system();

    int64_t mm_coeff_size(uint32_t id);
    int64_t mm_shape_size(uint32_t id);
    double mm_sigma(uint32_t id, int64_t idx);
    std::vector<float> mm_coeff(uint32_t id);
    std::vector<float> mm_basic_shape_coeff(uint32_t id, int64_t shape_id);
    void mm_coeff_set(uint32_t id, const std::vector<float>& coeff);

    int64_t mm_color_coeff_size(uint32_t id);
    int64_t mm_color_shape_size(uint32_t id);
    double mm_color_sigma(uint32_t id, int64_t idx);
    std::vector<float> mm_color_coeff(uint32_t id);
    std::vector<float> mm_color_basic_shape_coeff(uint32_t id, int64_t shape_id);
    void mm_color_coeff_set(uint32_t id, const std::vector<float>& coeff);

    std::vector<jtk::vec3<uint8_t>> mesh_texture_to_vertexcolors(uint32_t id);

    int64_t mm_to_mesh(int32_t id);

    void set_vertex_colors(uint32_t id, const std::vector<jtk::vec3<uint8_t>>& colors);

    bool vertices_to_csv(int64_t id, const char* filename);

    bool triangles_to_csv(int64_t id, const char* filename);

    void hide(int64_t id);

    void show(int64_t id);

    void hide();

    void show();

    void set_matcap(int64_t id, int64_t clr_id);

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

    uint32_t get_index(int x, int y);

    uint32_t get_id(int x, int y);

    int64_t marching_cubes(const jtk::boundingbox3d<float>& bb, uint64_t width, uint64_t height, uint64_t depth, float isovalue, double(*fun_ptr)(double, double, double));

    std::vector<jtk::vec3<uint32_t>> triangles(uint32_t id);

    std::vector<jtk::vec3<float>> vertices(uint32_t id);

    bool write(uint32_t id, const char* filename);

    void load_face_detector(const char* filename);

    void set_show_face_detector(bool b);

    std::vector<std::pair<long, long>> face_detector_predict();

    void fit_mm_to_mesh(uint32_t mm_id, uint32_t mesh_id);

    void fit_mm_to_partial_positions(uint32_t mm_id, const std::vector<uint32_t>& vertex_indices, const std::vector<jtk::vec3<float>>& vertex_positions);

    void fit_mm(uint32_t mm_id, uint32_t mesh_id, const std::vector<uint32_t>& vertex_indices, const std::vector<jtk::vec3<float>>& vertex_positions);
    
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

    void render_scene();

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
    std::unique_ptr<face_detector> p_face_detector;

    std::mutex _mut;

    bool _show_face_detector;
  };