#pragma once

#include <algorithm>
#include <unordered_map>
#include <vector>

template <typename _RandomAccessIterator, typename _Compare>
inline void _up_heap(_RandomAccessIterator _first, _RandomAccessIterator _pos, _Compare _comp)
  {
  auto _index = _pos - _first;
  auto _parent = (_index - 1) / 2;
  auto _val = *_pos;

  while (_index > 0 && _comp(*(_first + _parent), _val)) {
    *(_first + _index) = *(_first + _parent);
    _index = _parent;
    _parent = (_parent - 1) / 2;
    }

  if (_pos != (_first + _index))
    *(_first + _index) = _val;
  }

/////////////////////////////////////////////////////////////////////////////////

template <typename _RandomAccessIterator, typename _Compare>
inline void _down_heap(_RandomAccessIterator _first,
  _RandomAccessIterator _last,
  _RandomAccessIterator _pos,
  _Compare _comp)
  {
  auto _len = _last - _first;
  auto _index = _pos - _first;
  auto _left = _index * 2 + 1;
  auto _right = _index * 2 + 2;
  auto _largest = _right;
  auto _val = *_pos;

  while (_index < _len) {
    if (_right < _len) {
      _largest = _comp(*(_first + _right), *(_first + _left)) ? _left : _right;
      }
    else if (_left < _len) {
      _largest = _left;
      }
    else {
      // Force termination
      _largest = _len;
      }

    if (_largest < _len && _comp(_val, *(_first + _largest))) {
      *(_first + _index) = *(_first + _largest);
      _index = _largest;
      _left = _index * 2 + 1;
      _right = _index * 2 + 2;
      }
    else
      break;
    }

  if (_pos != (_first + _index))
    *(_first + _index) = _val;
  }

/////////////////////////////////////////////////////////////////////////////////

template <typename _RandomAccessIterator, typename _Compare>
inline void update_heap(_RandomAccessIterator _first,
  _RandomAccessIterator _last,
  _RandomAccessIterator _pos,
  _Compare _comp)
  {
  auto _index = (_pos - _first);
  auto _parent = (_index - 1) / 2;

  if (_index > 0 && _comp(*(_first + _parent), *(_pos)))
    _up_heap(_first, _pos, _comp);
  else
    _down_heap(_first, _last, _pos, _comp);
  }

/////////////////////////////////////////////////////////////////////////////////

template <class Key,
  class Data,
  class HashKey = std::hash<Key>,
  class KeyEquality = std::equal_to<Key>,
  class CompareData = std::less<Data>>
  class hashed_heap
  {
  private:
    using hash_type = std::unordered_map<Key, Data, HashKey, KeyEquality>;
    using heap_type = std::vector<Key>;

  public:
    using my_type = hashed_heap<Key, Data, HashKey, KeyEquality, CompareData>;
    using value_type = std::pair<Key, Data>;
    using const_reference = const std::pair<Key, Data>&;

  public:
    hashed_heap()
      : m_hash()
      , m_heap()
      {}

    /////////////////////////////////////////////////////////////////////////////////

    hashed_heap(const hashed_heap& i_other)
      : m_hash(i_other.m_hash)
      , m_heap(i_other.m_heap)
      {}

    /////////////////////////////////////////////////////////////////////////////////

    my_type& operator=(const my_type& i_right)
      {
      my_type temp(i_right);
      swap(temp);
      return (*this);
      }

    /////////////////////////////////////////////////////////////////////////////////

    void swap(my_type& i_right)
      {
      m_hash.swap(i_right.m_hash);
      m_heap.swap(i_right.m_heap);
      }

    /////////////////////////////////////////////////////////////////////////////////

    void reserve(std::size_t i_capacity) { m_heap.reserve(i_capacity); }

    /////////////////////////////////////////////////////////////////////////////////

    void push(const_reference i_val)
      {
      auto it = m_hash.find(i_val.first);
      if (it != m_hash.end() &&
        CompareData {}(i_val.second, it->second)) // element was already present, but with higher priority
        {
        it->second = i_val.second;
        auto iter =
          std::find_if(m_heap.begin(), m_heap.end(), [=](const Key& key) { return KeyEquality{}(key, i_val.first); });
        update_heap(m_heap.begin(), m_heap.end(), iter, Compare(m_hash, m_compare_data));
        }
      else if (it == m_hash.end()) {
        m_hash[i_val.first] = i_val.second;
        m_heap.push_back(i_val.first);
        update_heap(m_heap.begin(),
          m_heap.end(),
          m_heap.end() - 1,
          Compare(m_hash, m_compare_data)); // do not use std::push_heap, the STL does not guarantee that heaps
                                            // are implemented corresponding to update_heap
        }
      }

    /////////////////////////////////////////////////////////////////////////////////

    const hash_type& get_hash() const { return m_hash; }

    /////////////////////////////////////////////////////////////////////////////////

    const heap_type& get_heap() const { return m_heap; }

    /////////////////////////////////////////////////////////////////////////////////

    value_type top() const
      {
      std::pair<Key, Data> top_element(m_heap.front(), m_hash.at(m_heap.front()));
      return top_element;
      }

    /////////////////////////////////////////////////////////////////////////////////

    void pop() // do not use std::pop_heap, the STL does not guarantee that heaps are implemented corresponding to update_heap
      {
      Key k = m_heap.front();
      std::swap(m_heap.front(), m_heap.back());
      m_heap.pop_back();
      if (!m_heap.empty())
        update_heap(m_heap.begin(), m_heap.end(), m_heap.begin(), Compare(m_hash, m_compare_data));
      m_hash.erase(k);
      }

    /////////////////////////////////////////////////////////////////////////////////

    bool empty() const { return m_heap.empty(); }

    /////////////////////////////////////////////////////////////////////////////////

  private:
    class Compare
      {
      private:
        hash_type& m_hash;
        CompareData& m_compare_data;

      public:

        Compare(const Compare& other)
          : m_hash(other.m_hash)
          , m_compare_data(other.m_compare_data)
          {}

        void swap(Compare& other)
          {
          std::swap(m_hash, other.m_hash);
          std::swap(m_compare_data, other.m_compare_data);
          }

        Compare& operator=(const Compare& other)
          {
          Compare temp(other);
          swap(temp);
          return *this;
          }

        Compare(hash_type& i_hash, CompareData& i_compare_data)
          : m_hash(i_hash)
          , m_compare_data(i_compare_data)
          {}

        bool operator()(const Key& left, const Key& right) const { return m_compare_data(m_hash[right], m_hash[left]); }
      };

  private:
    hash_type m_hash;
    heap_type m_heap;
    CompareData m_compare_data;
  };