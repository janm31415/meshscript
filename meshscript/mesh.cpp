#include "mesh.h"
#include "io.h"

#include "primitives.h"

#include <jtk/geometry.h>

#include <algorithm>

#include <jtk/file_utils.h>

#include <trico/trico/trico.h>

extern "C"
  {
#include <libcork/triangle.h>
  }

#include <stb_image.h>

#include <iostream>
#include <list>

using namespace jtk;


jtk::image<uint32_t> make_dummy_texture(int w, int h, int block_size)
  {
  jtk::image<uint32_t> im(w, h);
  for (int y = 0; y < h; ++y)
    {
    uint32_t* p_clr = im.row(y);
    bool y_even = ((y / block_size) & 1) == 1;
    for (int x = 0; x < w; ++x, ++p_clr)
      {
      bool x_even = ((x / block_size) & 1) == 1;
      bool black = (x_even && y_even) || (!x_even && !y_even);
      *p_clr = black ? 0xff000000 : 0xffffffff;
      }
    }
  return im;
  }


void compute_bb(vec3<float>& min, vec3<float>& max, uint32_t nr_of_vertices, const vec3<float>* vertices)
  {
  if (nr_of_vertices == 0)
    return;
  min[0] = vertices[0][0];
  min[1] = vertices[0][1];
  min[2] = vertices[0][2];
  max[0] = min[0];
  max[1] = min[1];
  max[2] = min[2];
  for (uint32_t i = 1; i < nr_of_vertices; ++i)
    {
    min[0] = std::min<float>(min[0], vertices[i][0]);
    min[1] = std::min<float>(min[1], vertices[i][1]);
    min[2] = std::min<float>(min[2], vertices[i][2]);
    max[0] = std::max<float>(max[0], vertices[i][0]);
    max[1] = std::max<float>(max[1], vertices[i][1]);
    max[2] = std::max<float>(max[2], vertices[i][2]);
    }
  }

bool read_from_file(mesh& m, const std::string& filename)
  {
  std::string ext = jtk::get_extension(filename);
  if (ext.empty())
    return false;
  std::transform(ext.begin(), ext.end(), ext.begin(), [](char ch) {return (char)::tolower(ch); });
  if (ext == "stl")
    {
    if (!read_stl(m.vertices, m.triangles, filename.c_str()))
      return false;
    }
  else if (ext == "ply")
    {
    std::vector<jtk::vec3<float>> vertex_normals;
    std::vector<uint32_t> vertex_colors;
    if (!read_ply(filename.c_str(), m.vertices, vertex_normals, vertex_colors, m.triangles, m.uv_coordinates))
      return false;
    if (!vertex_colors.empty())
      {
      m.vertex_colors = convert_vertex_colors(vertex_colors);
      }
    }
  else if (ext == "off")
    {
    if (!read_off(m.vertices, m.triangles, filename.c_str()))
      return false;
    }
  else if (ext == "obj")
    {
    std::vector<jtk::vec3<float>> vertex_normals;
    std::vector<uint32_t> vertex_colors;
    if (!read_obj(filename.c_str(), m.vertices, vertex_normals, vertex_colors, m.triangles, m.uv_coordinates, m.texture))
      return false;
    if (!vertex_colors.empty())
      {
      m.vertex_colors = convert_vertex_colors(vertex_colors);
      }
    }
  else if (ext == "trc")
    {
    std::vector<jtk::vec3<float>> vertex_normals;
    std::vector<uint32_t> vertex_colors;
    if (!read_trc(filename.c_str(), m.vertices, vertex_normals, vertex_colors, m.triangles, m.uv_coordinates))
      return false;
    if (!vertex_colors.empty())
      {
      m.vertex_colors = convert_vertex_colors(vertex_colors);
      }
    }
  else
    return false;
  if (!m.uv_coordinates.empty() && (m.texture.width() == 0 || m.texture.height() == 0))
    m.texture = make_dummy_texture(512, 512);
  m.cs = get_identity();
  m.visible = true;
  return true;
  }

bool vertices_to_csv(const mesh& m, const std::string& filename)
  {
  std::vector<std::vector<std::string>> data;
  for (const auto& vertex : m.vertices)
    {
    std::vector<std::string> line;
    line.push_back(std::to_string(vertex[0]));
    line.push_back(std::to_string(vertex[1]));
    line.push_back(std::to_string(vertex[2]));
    data.push_back(line);
    }
  return csv_write(data, filename.c_str(), ",");
  }

bool triangles_to_csv(const mesh& m, const std::string& filename)
  {
  std::vector<std::vector<std::string>> data;
  for (const auto& tria : m.triangles)
    {
    std::vector<std::string> line;
    line.push_back(std::to_string(tria[0]));
    line.push_back(std::to_string(tria[1]));
    line.push_back(std::to_string(tria[2]));
    data.push_back(line);
    }
  return csv_write(data, filename.c_str(), ",");
  }

std::vector<uint32_t> convert_vertex_colors(const std::vector<jtk::vec3<float>>& vertex_colors)
  {
  std::vector<uint32_t> colors;
  colors.reserve(vertex_colors.size());
  for (auto& clr : vertex_colors)
    {
    uint32_t r = (uint32_t)(clr[0] * 255.f);
    uint32_t g = (uint32_t)(clr[1] * 255.f);
    uint32_t b = (uint32_t)(clr[2] * 255.f);
    r = r > 255 ? 255 : r;
    g = g > 255 ? 255 : g;
    b = b > 255 ? 255 : b;
    colors.push_back(0xff000000 | (b << 16) | (g << 8) | r);
    }
  return colors;
  }

std::vector<jtk::vec3<float>> convert_vertex_colors(const std::vector<uint32_t>& vertex_colors)
  {
  std::vector<jtk::vec3<float>> out;
  out.reserve(vertex_colors.size());
  for (uint32_t clr : vertex_colors)
    {
    uint32_t red = clr & 255;
    uint32_t green = (clr >> 8) & 255;
    uint32_t blue = (clr >> 16) & 255;
    out.emplace_back((float)red / 255.f, (float)green / 255.f, (float)blue / 255.f);
    }
  return out;
  }

bool write_to_file(const mesh& m, const std::string& filename)
  {
  std::string ext = jtk::get_extension(filename);
  if (ext.empty())
    return false;
  std::transform(ext.begin(), ext.end(), ext.begin(), [](char ch) {return (char)::tolower(ch); });
  if (ext == "stl")
    {
    return jtk::write_stl(m.vertices.data(), (uint32_t)m.triangles.size(), m.triangles.data(), nullptr, nullptr, filename.c_str());
    }
  else if (ext == "ply")
    {
    std::vector<uint32_t> colors = convert_vertex_colors(m.vertex_colors);
    std::vector<jtk::vec3<float>> normals;
    return write_ply(filename.c_str(), m.vertices, normals, colors, m.triangles, m.uv_coordinates);
    }
  else if (ext == "off")
    {
    return jtk::write_off((uint32_t)m.vertices.size(), m.vertices.data(), (uint32_t)m.triangles.size(), m.triangles.data(), filename.c_str());
    }
  else if (ext == "trc")
    {
    std::vector<uint32_t> colors = convert_vertex_colors(m.vertex_colors);
    std::vector<jtk::vec3<float>> normals;
    return write_trc(filename.c_str(), m.vertices, normals, colors, m.triangles, m.uv_coordinates);
    }
  else if (ext == "obj")
    {
    std::vector<uint32_t> colors = convert_vertex_colors(m.vertex_colors);
    std::vector<jtk::vec3<float>> normals;
    return write_obj(filename.c_str(), m.vertices, normals, colors, m.triangles, m.uv_coordinates, m.texture);
    }

  return false;
  }

void info(const mesh& m)
  {
  std::cout << "---------------------------------------" << std::endl;
  std::cout << "MESH" << std::endl;
  std::cout << "Triangles: " << m.triangles.size() << std::endl;
  std::cout << "Vertices: " << m.vertices.size() << std::endl;
  std::cout << "Coordinate system: " << std::endl;
  for (int i = 0; i < 4; ++i)
    {
    for (int j = 0; j < 4; ++j)
      {
      std::cout << m.cs[i + 4 * j] << " ";
      }
    std::cout << std::endl;
    }
  std::cout << "Vertex colors: " << (m.vertex_colors.empty() ? "No" : "Yes") << std::endl;
  std::cout << "UV coordinates: " << (m.uv_coordinates.empty() ? "No" : "Yes") << std::endl;
  std::cout << "Texture dimensions: " << m.texture.width() << " x " << m.texture.height() << std::endl;
  std::cout << "Visible: " << (m.visible ? "Yes" : "No") << std::endl;
  std::cout << "---------------------------------------" << std::endl;
  }

void cs_apply(mesh& m)
  {
  for (auto& v : m.vertices)
    {
    jtk::float4 V(v[0], v[1], v[2], 1.f);
    V = jtk::matrix_vector_multiply(m.cs, V);
    v[0] = V[0];
    v[1] = V[1];
    v[2] = V[2];
    }
  m.cs = jtk::get_identity();
  }


void make_cube(mesh& m, float w, float h, float d)
  {
  m = mesh();
  cube c;
  m.vertices = c.vertices;
  m.triangles = c.triangles;

  for (auto& v : m.vertices)
    {
    v[0] *= w;
    v[1] *= h;
    v[2] *= d;
    }

  m.cs = get_identity();
  m.visible = true;
  }

void make_sphere(mesh& m, float r, uint32_t subdivision_levels)
  {
  m = mesh();

  icosahedron ico;
  m.vertices = ico.vertices;
  m.triangles = ico.triangles;

  for (uint32_t iter = 0; iter < subdivision_levels; ++iter)
    jtk::dyadic_subdivide(m.vertices, m.triangles);

  for (auto& v : m.vertices)
    {
    v = (r*v) / jtk::length(v);
    }

  m.cs = get_identity();
  m.visible = true;
  }

void make_icosahedron(mesh& m, float r)
  {
  m = mesh();

  icosahedron ico;
  m.vertices = ico.vertices;
  m.triangles = ico.triangles;

  for (auto& v : m.vertices)
    {
    v = (r*v) / jtk::length(v);
    }

  m.cs = get_identity();
  m.visible = true;
  }

void make_cylinder(mesh& m, float r, float d, uint32_t n)
  {
  m = mesh();

  cylinder c(n);
  m.vertices = c.vertices;
  m.triangles = c.triangles;

  for (auto& v : m.vertices)
    {
    float rr = std::sqrt(v[0] * v[0] + v[1] * v[1]);
    if (rr)
      {
      v[0] *= r / rr;
      v[1] *= r / rr;
      }
    v[2] *= d;
    }

  m.cs = get_identity();
  m.visible = true;
  }

void triangulate_points(mesh& m, const std::vector<jtk::vec2<float>>& pts)
  {
  m = mesh();

  struct triangulateio in, out;

  in.numberofpoints = (int)pts.size();
  in.numberofpointattributes = 0;
  in.pointlist = new REAL[in.numberofpoints * 2];
  in.pointattributelist = nullptr;
  in.pointmarkerlist = new int[in.numberofpoints];
  for (int i = 0; i < in.numberofpoints; ++i)
    {
    in.pointlist[i * 2 + 0] = pts[i][0];
    in.pointlist[i * 2 + 1] = pts[i][1];
    in.pointmarkerlist[i] = 1;
    }
  in.numberofsegments = (int)pts.size();
  in.numberofholes = 0;
  in.numberofregions = 0;
  in.segmentlist = new int[in.numberofsegments * 2];
  in.segmentmarkerlist = new int[in.numberofsegments];
  for (int i = 0; i < in.numberofsegments; ++i)
    in.segmentmarkerlist[i] = 1;
  for (int i = 0; i < in.numberofsegments; ++i)
    {
    in.segmentlist[i * 2] = i;
    in.segmentlist[i * 2 + 1] = (i + 1) % pts.size();
    }
  in.numberoftriangles = 0;
  in.numberoftriangleattributes = 0;

  out.pointlist = nullptr;
  out.pointattributelist = nullptr;
  out.pointmarkerlist = nullptr;
  out.trianglelist = nullptr;
  out.segmentlist = nullptr;
  out.segmentmarkerlist = nullptr;

  char* params = (char*)("pzQYY");

  triangulate(params, &in, &out, nullptr);

  for (int i = 0; i < out.numberofpoints; ++i)
    {
    double x = out.pointlist[i * 2];
    double y = out.pointlist[i * 2 + 1];
    m.vertices.emplace_back((float)x, (float)y, 0.f);
    }
  for (int i = 0; i < out.numberoftriangles; ++i)
    {
    int a = out.trianglelist[i * 3];
    int b = out.trianglelist[i * 3 + 1];
    int c = out.trianglelist[i * 3 + 2];
    m.triangles.emplace_back((uint32_t)a, (uint32_t)b, (uint32_t)c);
    }


  m.cs = get_identity();
  m.visible = true;
  }

namespace
  {
  // https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order#:~:text=If%20the%20determinant%20is%20negative,q%20and%20r%20are%20collinear.
  bool polyline_is_clockwise(const std::vector<jtk::vec2<float>>& pts)
    {
    float sum = 0.f;
    for (size_t i = 0; i < pts.size(); ++i)
      {
      auto a = pts[i];
      auto b = pts[(i + 1) % pts.size()];
      float term = (b[0] - a[0])*(b[1] + a[1]);
      sum += term;
      }
    return sum > 0.f;
    }

  }

void revolve_points(mesh& m, const std::vector<jtk::vec2<float>>& pts, uint32_t n, bool closed)
  {
  bool clockwise = polyline_is_clockwise(pts);

  m = mesh();
  const double twopi = 2.0* 3.14159265358979323846264338327950288419716939937510;

  m.vertices.reserve(pts.size() * n);

  for (uint32_t i = 0; i < (uint32_t)pts.size(); ++i)
    {
    uint32_t idx = clockwise ? i : (uint32_t)pts.size() - 1 - i;
    for (uint32_t j = 0; j < n; ++j)
      {      
      const double angle = (double)j * twopi / (double)n ;
      m.vertices.emplace_back(pts[idx][0], (float)(pts[idx][1] * std::sin(angle)), (float)(pts[idx][1] * std::cos(angle)));
      }
    }

  uint32_t sz = (uint32_t)m.vertices.size();

  uint32_t top = closed ? (uint32_t)pts.size() : (uint32_t)pts.size() - 1;

  for (uint32_t i = 0; i < top; ++i)
    {
    for (uint32_t j = 0; j < n; ++j)
      {
      m.triangles.emplace_back(i*n + j, (((i + 1) % pts.size())*n + j), (i*n + (j + 1) % n));
      m.triangles.emplace_back((i*n + (j + 1) % n), (((i + 1) % pts.size())*n + j), (((i + 1) % pts.size())*n + (j + 1) % n));
      }
    }


  m.cs = get_identity();
  m.visible = true;
  }

void extrude_points(mesh& m, const std::vector<jtk::vec2<float>>& pts, float h)
  {
  m = mesh();

  bool clockwise = polyline_is_clockwise(pts);

  struct triangulateio in, out;

  in.numberofpoints = (int)pts.size();
  in.numberofpointattributes = 0;
  in.pointlist = new REAL[in.numberofpoints * 2];
  in.pointattributelist = nullptr;
  in.pointmarkerlist = new int[in.numberofpoints];
  for (int i = 0; i < in.numberofpoints; ++i)
    {
    in.pointlist[i * 2 + 0] = clockwise ? pts[i][0] : pts[in.numberofpoints - 1 - i][0];
    in.pointlist[i * 2 + 1] = clockwise ? pts[i][1] : pts[in.numberofpoints - 1 - i][1];
    in.pointmarkerlist[i] = 1;
    }
  in.numberofsegments = (int)pts.size();
  in.numberofholes = 0;
  in.numberofregions = 0;
  in.segmentlist = new int[in.numberofsegments * 2];
  in.segmentmarkerlist = new int[in.numberofsegments];
  for (int i = 0; i < in.numberofsegments; ++i)
    in.segmentmarkerlist[i] = 1;
  for (int i = 0; i < in.numberofsegments; ++i)
    {
    in.segmentlist[i * 2] = i;
    in.segmentlist[i * 2 + 1] = (i + 1) % pts.size();
    }
  in.numberoftriangles = 0;
  in.numberoftriangleattributes = 0;

  out.pointlist = nullptr;
  out.pointattributelist = nullptr;
  out.pointmarkerlist = nullptr;
  out.trianglelist = nullptr;
  out.segmentlist = nullptr;
  out.segmentmarkerlist = nullptr;

  char* params = (char*)("pzQYY");

  triangulate(params, &in, &out, nullptr);

  for (int i = 0; i < out.numberofpoints; ++i)
    {
    double x = out.pointlist[i * 2];
    double y = out.pointlist[i * 2 + 1];
    m.vertices.emplace_back((float)x, (float)y, 0.f);
    m.vertices.emplace_back((float)x, (float)y, h);
    }
  for (int i = 0; i < out.numberoftriangles; ++i)
    {
    int a = out.trianglelist[i * 3];
    int b = out.trianglelist[i * 3 + 1];
    int c = out.trianglelist[i * 3 + 2];
    m.triangles.emplace_back((uint32_t)c * 2, (uint32_t)b * 2, (uint32_t)a * 2);
    m.triangles.emplace_back((uint32_t)a * 2 + 1, (uint32_t)b * 2 + 1, (uint32_t)c * 2 + 1);
    }

  for (int i = 0; i < out.numberofsegments; ++i)
    {
    int a = out.segmentlist[2 * i];
    int b = out.segmentlist[2 * i + 1];
    m.triangles.emplace_back((uint32_t)b * 2, (uint32_t)a * 2, (uint32_t)a * 2 + 1);
    m.triangles.emplace_back((uint32_t)b * 2, (uint32_t)a * 2 + 1, (uint32_t)b * 2 + 1);
    }

  m.cs = get_identity();
  m.visible = true;
  }

inline std::vector<std::vector<uint32_t>> ordered_one_ring_vertices_from_vertex(uint32_t vertex_index, const adjacency_list& adj_list, const vec3<uint32_t>* triangles, bool oriented = true)
    {
    auto it = adj_list.begin(vertex_index);
    const auto it_end = adj_list.end(vertex_index);
    if (it == it_end)
      return std::vector<std::vector<uint32_t>>();
      
    std::list<uint32_t> trias;
    for (auto it2 = it; it2 != it_end; ++it2)
      {
      trias.push_back(*it2);
      }
    
    std::vector<std::vector<uint32_t>> batches;
    
    while (!trias.empty())
      {
      std::vector<uint32_t> current_loop;
      const auto tria = triangles[trias.front()];
      trias.pop_front();
      uint32_t v = 0;
      if (tria[v] != vertex_index)
        {
        ++v;
        if (tria[v] != vertex_index)
          ++v;
        }
      current_loop.push_back(tria[(v + 1) % 3]);
      current_loop.push_back(tria[(v + 2) % 3]);
      bool done = false;
      while (!done)
        {
        done = true;
        for (auto it3 = trias.begin(); it3 != trias.end(); ++it3)
          {
          const auto local_tria = triangles[*it3];
          v = 0;
          if (local_tria[v] != vertex_index)
            {
            ++v;
            if (local_tria[v] != vertex_index)
              ++v;
            }
          if (local_tria[(v + 1) % 3] == current_loop.back())
            {
            done = false;
            current_loop.push_back(local_tria[(v + 2) % 3]);
            trias.erase(it3);
            break;
            }
          else if (local_tria[(v + 2) % 3] == current_loop.front())
            {
            done = false;
            current_loop.insert(current_loop.begin(), local_tria[(v + 1) % 3]);
            trias.erase(it3);
            break;
            }
          else if (!oriented && (local_tria[(v + 2) % 3] == current_loop.back()))
            {
            done = false;
            current_loop.push_back(local_tria[(v + 1) % 3]);
            trias.erase(it3);
            break;
            }
          else if (!oriented && (local_tria[(v + 1) % 3] == current_loop.front()))
            {
            done = false;
            current_loop.insert(current_loop.begin(), local_tria[(v + 2) % 3]);
            trias.erase(it3);
            break;
            }
            
          }
        }
      if (current_loop.back() == current_loop.front())
        current_loop.pop_back();
      batches.push_back(current_loop);
      }
    return batches;
    }
    
void butterfly(std::vector<jtk::vec3<float>>& vertices, std::vector<jtk::vec3<uint32_t>>& triangles)
  {
  using namespace jtk;
  
  const uint32_t nr_of_vertices = vertices.size();
  const uint32_t nr_of_triangles = triangles.size();

  assert(nr_of_triangles % 4 == 0);

  uint32_t original_nr_of_vertices = triangles.front()[0];
  const uint32_t original_nr_of_triangles = nr_of_triangles / 4;
  for (uint32_t t = 0; t < original_nr_of_triangles; ++t)
    {
    original_nr_of_vertices = std::min<uint32_t>(original_nr_of_vertices, triangles[t][0]);
    original_nr_of_vertices = std::min<uint32_t>(original_nr_of_vertices, triangles[t][1]);
    original_nr_of_vertices = std::min<uint32_t>(original_nr_of_vertices, triangles[t][2]);
    }
    
  std::vector<jtk::vec3<uint32_t>> old_triangles;
  old_triangles.reserve(original_nr_of_triangles);
  for (uint32_t t = 0; t < original_nr_of_triangles; ++t)
    {
    const uint32_t t0 = original_nr_of_triangles + t * 3;
    const uint32_t t1 = original_nr_of_triangles + t * 3 + 1;
    const uint32_t t2 = original_nr_of_triangles + t * 3 + 2;
    old_triangles.emplace_back(triangles[t0][0], triangles[t1][1], triangles[t2][2]);
    }
    
  adjacency_list new_adj_list(nr_of_vertices, triangles.data(), nr_of_triangles);
  adjacency_list old_adj_list(original_nr_of_vertices, old_triangles.data(), original_nr_of_triangles);
  for (uint32_t v = original_nr_of_vertices; v < (uint32_t)vertices.size(); ++v)
    {
    auto neighbours = one_ring_vertices_from_vertex(v, new_adj_list, triangles.data());
    uint32_t e[2];
    int e_index = 0;
    for (const auto& n : neighbours)
      {
      if (n < original_nr_of_vertices)
        e[e_index++] = n;
      }
      
    if (is_boundary_edge(e[0], e[1], old_adj_list))
      {
      uint32_t e2[2];
      for (int k = 0; k < 2; ++k)
        {
        auto neighbours2 = one_ring_vertices_from_vertex(e[k], old_adj_list, old_triangles.data());
        for (const auto& n : neighbours2)
          {
          if (n != e[(k+1)%2] && is_boundary_edge(e[k], n, old_adj_list))
            {
            e2[k] = n;
            break;
            }
          }
        }
      auto new_vertex_position = (vertices[e2[0]] + vertices[e2[1]])/-16.f + 9.f*(vertices[e[0]] + vertices[e[1]])/16.f;
      vertices[v] = new_vertex_position;
      }
    else
      {
      //auto valence0 = one_ring_vertices_from_vertex(e[0], old_adj_list, old_triangles.data()).size();
      //auto valence1 = one_ring_vertices_from_vertex(e[1], old_adj_list, old_triangles.data()).size();
      auto n0 = ordered_one_ring_vertices_from_vertex(e[0], old_adj_list, old_triangles.data()).front();
      auto n1 = ordered_one_ring_vertices_from_vertex(e[1], old_adj_list, old_triangles.data()).front();
      
      size_t e1_position = std::distance(n0.begin(), std::find(n0.begin(), n0.end(), e[1]));
      size_t e0_position = std::distance(n1.begin(), std::find(n1.begin(), n1.end(), e[0]));
      
      assert(e1_position < n0.size());
      assert(e0_position < n1.size());
      
      if (n0.size() == 6 && n1.size() == 6) // regular butterfly
        {
        uint32_t b0_up = n0[(e1_position+1)%n0.size()];
        uint32_t b0_down = n0[(e1_position+n0.size()-1)%n0.size()];
        uint32_t c0_up = n0[(e1_position+2)%n0.size()];
        uint32_t c0_down = n0[(e1_position+n0.size()-2)%n0.size()];
        
        uint32_t b1_down = n1[(e0_position+1)%n1.size()];
        uint32_t b1_up = n1[(e0_position+n1.size()-1)%n1.size()];
        uint32_t c1_down = n1[(e0_position+2)%n1.size()];
        uint32_t c1_up = n1[(e0_position+n0.size()-2)%n1.size()];
        
        assert(b0_up == b1_up);
        assert(b0_down == b1_down);
        
        auto new_vertex_position = (vertices[e[0]] + vertices[e[1]])*0.5f + (vertices[b0_up] + vertices[b0_down])/8.f - (vertices[c0_up] + vertices[c0_down] + vertices[c1_up] + vertices[c1_down])/16.f;
        vertices[v] = new_vertex_position;
        }
      else if (n0.size() == 6 || n1.size() == 6)
        {
        if (n0.size() == 6)
          {
          std::swap(n0, n1);
          std::swap(e0_position, e1_position);
          }
        std::vector<float> weights;
        weights.reserve(n0.size());
        if (n0.size() == 3)
          {
          weights.push_back(5.f/12.f);
          weights.push_back(-1.f/12.f);
          weights.push_back(-1.f/12.f);
          }
        else if (n0.size() == 4)
          {
          weights.push_back(3.f/8.f);
          weights.push_back(0.f);
          weights.push_back(-1.f/8.f);
          weights.push_back(0.f);
          }
        else
          {
          for (int j = 0; j < (int)n0.size(); ++j)
            {
            weights.push_back((0.25+std::cos(6.28318530718f*(float)j/(float)n0.size()) + 0.5*std::cos(12.5663706144f*(float)j/(float)n0.size()) )/(float)n0.size());
            }
          }
          
        float q = 1.f - std::accumulate(weights.begin(), weights.end(), 0.f);
        jtk::vec3<float> new_vertex_position = q*vertices[n1[e0_position]];
        for (int j = 0; j < (int)n0.size(); ++j)
          {
          new_vertex_position = new_vertex_position + vertices[n0[(j+e1_position)%n0.size()]]*weights[j];
          }
        vertices[v] = new_vertex_position;
        }
      else
        {
        std::vector<float> weights0;
        weights0.reserve(n0.size());
        if (n0.size() == 3)
          {
          weights0.push_back(5.f/12.f);
          weights0.push_back(-1.f/12.f);
          weights0.push_back(-1.f/12.f);
          }
        else if (n0.size() == 4)
          {
          weights0.push_back(3.f/8.f);
          weights0.push_back(0.f);
          weights0.push_back(-1.f/8.f);
          weights0.push_back(0.f);
          }
        else
          {
          for (int j = 0; j < (int)n0.size(); ++j)
            {
            weights0.push_back((0.25+std::cos(6.28318530718f*(float)j/(float)n0.size()) + 0.5*std::cos(12.5663706144f*(float)j/(float)n0.size()) )/(float)n0.size());
            }
          }
          
        float q0 = 1.f - std::accumulate(weights0.begin(), weights0.end(), 0.f);
        jtk::vec3<float> new_vertex_position_0 = q0*vertices[n1[e0_position]];
        for (int j = 0; j < (int)n0.size(); ++j)
          {
          new_vertex_position_0 = new_vertex_position_0 + vertices[n0[(j+e1_position)%n0.size()]]*weights0[j];
          }
          
        std::vector<float> weights1;
        weights1.reserve(n1.size());
        if (n1.size() == 3)
          {
          weights1.push_back(5.f/12.f);
          weights1.push_back(-1.f/12.f);
          weights1.push_back(-1.f/12.f);
          }
        else if (n1.size() == 4)
          {
          weights1.push_back(3.f/8.f);
          weights1.push_back(0.f);
          weights1.push_back(-1.f/8.f);
          weights1.push_back(0.f);
          }
        else
          {
          for (int j = 0; j < (int)n1.size(); ++j)
            {
            weights1.push_back((0.25+std::cos(6.28318530718f*(float)j/(float)n1.size()) + 0.5*std::cos(12.5663706144f*(float)j/(float)n1.size()) )/(float)n1.size());
            }
          }
          
        float q1 = 1.f - std::accumulate(weights1.begin(), weights1.end(), 0.f);
        jtk::vec3<float> new_vertex_position_1 = q1*vertices[n0[e1_position]];
        for (int j = 0; j < (int)n1.size(); ++j)
          {
          new_vertex_position_1 = new_vertex_position_1 + vertices[n1[(j+e0_position)%n1.size()]]*weights1[j];
          }
          
        vertices[v] = (new_vertex_position_0 + new_vertex_position_1)*0.5f;
        }
      }
    }
  }
