#include "db.h"
#include "mesh.h"
#include "mm.h"
#include "pc.h"
#include "sp.h"
#include "shape_predictor.h"
#include "im.h"

#include <cassert>

using namespace jtk;

db::db()
  {

  }

db::~db()
  {
  clear();
  }

void db::swap(db& other)
  {
  std::swap(meshes, other.meshes);
  std::swap(meshes_deleted, other.meshes_deleted);
  std::swap(pcs, other.pcs);
  std::swap(pcs_deleted, other.pcs_deleted);
  std::swap(mms, other.mms);
  std::swap(mms_deleted, other.mms_deleted);
  std::swap(sps, other.sps);
  std::swap(sps_deleted, other.sps_deleted);
  std::swap(images, other.images);
  std::swap(images_deleted, other.images_deleted);
  }

void db::create_mesh(mesh*& new_mesh, uint32_t& id)
  {
  mesh* m = new mesh();
  id = make_db_id(MESH_KEY, (uint32_t)meshes.size());
  new_mesh = m;
  meshes.push_back(std::make_pair(id, m));
  meshes_deleted.push_back(std::make_pair(id, nullptr));
  }

mesh* db::get_mesh(uint32_t id) const
  {
  auto key = get_db_key(id);
  auto vector_index = get_db_vector_index(id);
  if (key != MESH_KEY)
    return nullptr;
  if (vector_index >= meshes.size())
    return nullptr;
  assert(meshes[vector_index].first == id);
  return meshes[vector_index].second;
  }

bool db::is_mesh(uint32_t id) const
  {
  return get_db_key(id) == MESH_KEY;
  }

void db::create_pc(pc*& new_pointcloud, uint32_t& id)
  {
  pc* m = new pc();
  id = make_db_id(PC_KEY, (uint32_t)pcs.size());
  new_pointcloud = m;
  pcs.push_back(std::make_pair(id, m));
  pcs_deleted.push_back(std::make_pair(id, nullptr));
  }

pc* db::get_pc(uint32_t id) const
  {
  auto key = get_db_key(id);
  auto vector_index = get_db_vector_index(id);
  if (key != PC_KEY)
    return nullptr;
  if (vector_index >= pcs.size())
    return nullptr;
  assert(pcs[vector_index].first == id);
  return pcs[vector_index].second;
  }

bool db::is_pc(uint32_t id) const
  {
  return get_db_key(id) == PC_KEY;
  }

void db::create_mm(mm*& new_mm, uint32_t& id)
  {
  mm* m = new mm();
  id = make_db_id(MM_KEY, (uint32_t)mms.size());
  new_mm = m;
  mms.push_back(std::make_pair(id, m));
  mms_deleted.push_back(std::make_pair(id, nullptr));
  }

mm* db::get_mm(uint32_t id) const
  {
  auto key = get_db_key(id);
  auto vector_index = get_db_vector_index(id);
  if (key != MM_KEY)
    return nullptr;
  if (vector_index >= mms.size())
    return nullptr;
  assert(mms[vector_index].first == id);
  return mms[vector_index].second;
  }

bool db::is_mm(uint32_t id) const
  {
  return get_db_key(id) == MM_KEY;
  }


void db::create_sp(sp*& new_sp, uint32_t& id)
  {
  sp* m = new sp();
  id = make_db_id(SP_KEY, (uint32_t)sps.size());
  new_sp = m;
  sps.push_back(std::make_pair(id, m));
  sps_deleted.push_back(std::make_pair(id, nullptr));
  }

sp* db::get_sp(uint32_t id) const
  {
  auto key = get_db_key(id);
  auto vector_index = get_db_vector_index(id);
  if (key != SP_KEY)
    return nullptr;
  if (vector_index >= sps.size())
    return nullptr;
  assert(sps[vector_index].first == id);
  return sps[vector_index].second;
  }

bool db::is_sp(uint32_t id) const
  {
  return get_db_key(id) == SP_KEY;
  }


void db::create_image(im*& new_im, uint32_t& id)
  {
  im* m = new im();
  id = make_db_id(IM_KEY, (uint32_t)images.size());
  new_im = m;
  images.push_back(std::make_pair(id, m));
  images_deleted.push_back(std::make_pair(id, nullptr));
  }

im* db::get_image(uint32_t id) const
  {
  auto key = get_db_key(id);
  auto vector_index = get_db_vector_index(id);
  if (key != IM_KEY)
    return nullptr;
  if (vector_index >= images.size())
    return nullptr;
  assert(images[vector_index].first == id);
  return images[vector_index].second;
  }

bool db::is_image(uint32_t id) const
  {
  return get_db_key(id) == IM_KEY;
  }

void db::delete_object(uint32_t id)
  {
  auto key = get_db_key(id);
  auto vector_index = get_db_vector_index(id);
  switch (key)
    {
    case MESH_KEY:
      if (meshes[vector_index].second)
        {
        meshes_deleted[vector_index].second = meshes[vector_index].second;
        meshes[vector_index].second = nullptr;
        }          
      break;
    case PC_KEY:
      if (pcs[vector_index].second)
        {
        pcs_deleted[vector_index].second = pcs[vector_index].second;
        pcs[vector_index].second = nullptr;
        }
      break;
    case MM_KEY:
      if (mms[vector_index].second)
        {
        mms_deleted[vector_index].second = mms[vector_index].second;
        mms[vector_index].second = nullptr;
        }
      break;
    case SP_KEY:
      if (sps[vector_index].second)
        {
        sps_deleted[vector_index].second = sps[vector_index].second;
        sps[vector_index].second = nullptr;
        }
      break;
    case IM_KEY:
      if (images[vector_index].second)
        {
        images_deleted[vector_index].second = images[vector_index].second;
        images[vector_index].second = nullptr;
        }
      break;
    }
  }

void db::restore_object(uint32_t id)
  {
  auto key = get_db_key(id);
  auto vector_index = get_db_vector_index(id);
  switch (key)
    {
    case MESH_KEY:
      if (!meshes[vector_index].second)
        {
        meshes[vector_index].second = meshes_deleted[vector_index].second;
        meshes_deleted[vector_index].second = nullptr;
        }
      break;    
    case PC_KEY:
      if (!pcs[vector_index].second)
        {
        pcs[vector_index].second = pcs_deleted[vector_index].second;
        pcs_deleted[vector_index].second = nullptr;
        }
      break;
    case MM_KEY:
      if (!mms[vector_index].second)
        {
        mms[vector_index].second = mms_deleted[vector_index].second;
        mms_deleted[vector_index].second = nullptr;
        }
      break;
    case SP_KEY:
      if (!sps[vector_index].second)
        {
        sps[vector_index].second = sps_deleted[vector_index].second;
        sps_deleted[vector_index].second = nullptr;
        }
      break;
    case IM_KEY:
      if (!images[vector_index].second)
        {
        images[vector_index].second = images_deleted[vector_index].second;
        images_deleted[vector_index].second = nullptr;
        }
      break;
    }
  }

namespace
  {
  template <class T>
  void delete_objects(std::vector<std::pair<uint32_t, T*>>& vec)
    {
    for (auto& p_obj : vec)
      {
      delete p_obj.second;
      p_obj.second = nullptr;
      }
    vec.clear();
    }
  }

void db::clear()
  {
  delete_objects(meshes);
  delete_objects(meshes_deleted);
  delete_objects(pcs);
  delete_objects(pcs_deleted);
  delete_objects(mms);
  delete_objects(mms_deleted);
  delete_objects(sps);
  delete_objects(sps_deleted);
  delete_objects(images);
  delete_objects(images_deleted);
  }

std::vector<vec3<float>>* get_vertices(const db& _db, uint32_t id)
  {
  auto key = get_db_key(id);
  switch (key)
    {
    case MESH_KEY:
      return &_db.get_mesh(id)->vertices;
      break;
    case PC_KEY:
      return &_db.get_pc(id)->vertices;
      break;
    case MM_KEY:
      return &_db.get_mm(id)->vertices;
      break;
    }
  return nullptr;
  }

std::vector<vec3<uint32_t>>* get_triangles(const db& _db, uint32_t id)
  {
  auto key = get_db_key(id);
  switch (key)
    {
    case MESH_KEY:
      return &_db.get_mesh(id)->triangles;
      break;
    case MM_KEY:
      return &_db.get_mm(id)->shape.triangles;
      break;
    }
  return nullptr;
  }

float4x4* get_cs(const db& _db, uint32_t id)
  {
  auto key = get_db_key(id);
  switch (key)
    {
    case MESH_KEY:
      return &_db.get_mesh(id)->cs;
      break;
    case PC_KEY:
      return &_db.get_pc(id)->cs;
      break;
    case MM_KEY:
      return &_db.get_mm(id)->cs;
      break;
    }
  return nullptr;
  }

bool is_visible(const db& _db, uint32_t id)
  {
  auto key = get_db_key(id);
  switch (key)
    {
    case MESH_KEY:
      return _db.get_mesh(id)->visible;
      break;
    case PC_KEY:
      return _db.get_pc(id)->visible;
      break;
    case MM_KEY:
      return _db.get_mm(id)->visible;
      break;
    }
  return false;
  }