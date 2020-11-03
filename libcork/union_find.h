/*
https://en.wikipedia.org/wiki/Disjoint-set_data_structure

In computer science, a union-find data structure, also called a disjoint-set data structure or merge–find set,
is a data structure that keeps track of a set of elements partitioned into a number of disjoint (nonoverlapping)
subsets. It supports two useful operations:

Find: Determine which subset a particular element is in. Find typically returns an item from this set that serves
      as its "representative"; by comparing the result of two Find operations, one can determine whether two elements
      are in the same subset.
Union: Join two subsets into a single subset.

*/
#pragma once

#include <vector>
#include <stdint.h>

class union_find
  {
  public:
    union_find(uint32_t N) : ids(N), rank(N, 0)
      {
      for (uint32_t i = 0; i < N; ++i)
        ids[i] = i;
      }

    uint32_t find(uint32_t i)
      {
      uint32_t id = i;
      while (ids[id] != id)
        id = ids[id];
      ids[i] = id;
      return id;
      }

    uint32_t union_sets(uint32_t i, uint32_t j)
      {
      uint32_t iid = find(i);
      uint32_t jid = find(j);
      if (iid == jid)
        return iid;
      if (rank[iid] > rank[jid])
        {
        return ids[jid] = iid;
        }
      else if (rank[iid] < rank[jid])
        {
        return ids[iid] = jid;
        }
      else
        {
        ++rank[jid];
        return ids[jid] = iid;
        }
      }

  private:
    std::vector<uint32_t> ids;
    std::vector<uint32_t> rank;
  };

