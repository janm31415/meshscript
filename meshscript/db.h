#pragma once

#include <stdint.h>
#include <vector>

#include <jtk/vec.h>
#include <jtk/qbvh.h>

struct mesh;
struct pc;
struct mm;
struct sp;
struct im;
struct deform_tool;

#define MESH_KEY 1 // mesh
#define PC_KEY 2   // point cloud
#define MM_KEY 3   // morphable model
#define SP_KEY 4   // shape predictor
#define IM_KEY 5   // image
#define DT_KEY 6   // deform tool

inline uint32_t get_db_key(uint32_t id)
  {
  return (id >> 29);
  }

inline uint32_t get_db_vector_index(uint32_t id)
  {
  return (id & 0x1fffffff);
  }

inline uint32_t make_db_id(uint32_t key, uint32_t vector_index)
  {
  return (key << 29) | vector_index;
  }

class db
  {
  public:
    db();
    ~db();

    db(db const&) = delete;
    db(db&&) = delete;
    db& operator=(db const&) = delete;
    db& operator=(db&&) = delete;

    void swap(db& other);

    void create_mesh(mesh*& mesh, uint32_t& id);
    mesh* get_mesh(uint32_t id) const;
    bool is_mesh(uint32_t id) const;   

    void create_pc(pc*& pointcloud, uint32_t& id);
    pc* get_pc(uint32_t id) const;
    bool is_pc(uint32_t id) const;

    void create_mm(mm*& morph, uint32_t& id);
    mm* get_mm(uint32_t id) const;
    bool is_mm(uint32_t id) const;

    void create_sp(sp*& shape_pred, uint32_t& id);
    sp* get_sp(uint32_t id) const;
    bool is_sp(uint32_t id) const;

    void create_image(im*& i, uint32_t& id);
    im* get_image(uint32_t id) const;
    bool is_image(uint32_t id) const;

    void create_deform_tool(deform_tool*& i, uint32_t& id);
    deform_tool* get_deform_tool(uint32_t id) const;
    bool is_deform_tool(uint32_t id) const;

    void delete_object(uint32_t id);
    void restore_object(uint32_t id);
    
    void clear();

    const std::vector<std::pair<uint32_t, mesh*>>& get_meshes() const { return meshes; }
    const std::vector<std::pair<uint32_t, pc*>>& get_pcs() const { return pcs; }
    const std::vector<std::pair<uint32_t, mm*>>& get_mms() const { return mms; }
    const std::vector<std::pair<uint32_t, sp*>>& get_sps() const { return sps; }
    const std::vector<std::pair<uint32_t, im*>>& get_images() const { return images; }
    const std::vector<std::pair<uint32_t, deform_tool*>>& get_deform_tools() const { return deform_tools; }

  private:
    std::vector<std::pair<uint32_t, mesh*>> meshes, meshes_deleted;
    std::vector<std::pair<uint32_t, pc*>> pcs, pcs_deleted;
    std::vector<std::pair<uint32_t, mm*>> mms, mms_deleted;
    std::vector<std::pair<uint32_t, sp*>> sps, sps_deleted;
    std::vector<std::pair<uint32_t, im*>> images, images_deleted;
    std::vector<std::pair<uint32_t, deform_tool*>> deform_tools, deform_tools_deleted;
  };

std::vector<jtk::vec3<float>>* get_vertices(const db& _db, uint32_t id);
std::vector<jtk::vec3<uint32_t>>* get_triangles(const db& _db, uint32_t id);
jtk::float4x4* get_cs(const db& _db, uint32_t id);
bool is_visible(const db& _db, uint32_t id);