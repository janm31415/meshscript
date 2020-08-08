#pragma once

#include <jtk/image.h>

#include <stdint.h>
#include <map>
#include <vector>

struct matcap
  {
  jtk::image<uint32_t> im;
  uint32_t cavity_clr;
  };

void make_matcap_gray(matcap& _matcap);
void make_matcap_brown(matcap& _matcap);
void make_matcap_red_wax(matcap& _matcap);
void make_matcap_sketch(matcap& _matcap);

struct matcapmap
  {
  matcapmap();
  std::vector<matcap> matcaps;
  std::map<uint32_t, uint32_t> map_db_id_to_matcap;

  const matcap& get_matcap(uint32_t id) const
    {
    auto it = map_db_id_to_matcap.find(id);
    if (it == map_db_id_to_matcap.end())
      return matcaps.front();
    return matcaps[it->second];
    }
  };