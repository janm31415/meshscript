#pragma once

#include <jtk/qbvh.h>
#include <jtk/image.h>
#include <jtk/vec.h>

#include <stdint.h>

#include <memory>
#include <string>

struct mesh
  {
  std::vector<jtk::vec3<float>> vertices;
  std::vector<jtk::vec3<uint32_t>> triangles;
  std::vector<jtk::vec3<float>> vertex_colors;
  std::vector<jtk::vec3<jtk::vec2<float>>> uv_coordinates;
  jtk::image<uint32_t> texture;
  jtk::float4x4 cs;
  bool visible;
  };

jtk::image<uint32_t> make_dummy_texture(int w, int h, int block_size = 32);
std::vector<uint32_t> convert_vertex_colors(const std::vector<jtk::vec3<float>>& vertex_colors);
std::vector<jtk::vec3<float>> convert_vertex_colors(const std::vector<uint32_t>& vertex_colors);

void compute_bb(jtk::vec3<float>& min, jtk::vec3<float>& max, uint32_t nr_of_vertices, const jtk::vec3<float>* vertices);
bool read_from_file(mesh& m, const std::string& filename);

bool vertices_to_csv(const mesh& m, const std::string& filename);
bool triangles_to_csv(const mesh& m, const std::string& filename);

bool write_to_file(const mesh& m, const std::string& filename);

void info(const mesh& m);

void cs_apply(mesh& m);


void make_cube(mesh& m, float w, float h, float d);
void make_sphere(mesh& m, float r, uint32_t subdivision_levels);
void make_icosahedron(mesh& m, float r);
void make_cylinder(mesh& m, float r, float d, uint32_t n);

void triangulate_points(mesh& m, const std::vector<jtk::vec2<float>>& pts);
void extrude_points(mesh& m, const std::vector<jtk::vec2<float>>& pts, float h);
void revolve_points(mesh& m, const std::vector<jtk::vec2<float>>& pts, uint32_t n, bool closed);
