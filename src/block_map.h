
#ifndef __BLOCK_MAP__
#define __BLOCK_MAP__
#include <algorithm>
#include <utility>

template<typename Key, typename T, typename KeyUnmapFun>
struct block_map_iterator {
  
  typedef std::forward_iterator_tag iterator_category;
  typedef std::pair<Key, T> value_type;
  typedef ptrdiff_t difference_type;
  typedef size_t size_type;
  typedef const std::pair<Key, T>* pointer;
  typedef const std::pair<Key, T>& reference;
  typedef block_map_iterator<Key, T, KeyUnmapFun> iterator;
  

  
  KeyUnmapFun key_unmap;
    T* data;
    size_type pos;
    size_type len;
    T empty;
    value_type return_val;
  

  
  block_map_iterator() {}
    
  block_map_iterator(T* d, size_type p, size_type l, T e) : 
    data(d), pos(p), len(l), empty(e)
  {
    if (pos < len && data[pos] == empty) incr();
  }
    
  block_map_iterator(const iterator& x) : 
    data(x.data), pos(x.pos), len(x.len), empty(x.empty) { }
  
  iterator& operator=(const iterator& x) {
    data = x.data;
    pos = x.pos;
    len = x.len;
    empty = x.empty;
    return *this;
  }
  

  
  void incr() {
      while (pos < len && data[++pos] == empty) { }
  }
    
  block_map_iterator& operator++() { incr(); return *this; }
  block_map_iterator operator++(int) { 
    block_map_iterator tmp = *this; incr(); return tmp;
  }
  

  
  void update() {
    return_val.first = key_unmap(pos);
    return_val.second = data[pos];
  }
    
  reference operator*() {
    update();
    return return_val;
  }
  
  pointer operator->() {
    update();
    return &return_val;
  }
  

  
  bool operator==(const block_map_iterator& x) const {
    return data == x.data && pos == x.pos;
  }
  bool operator!=(const block_map_iterator& x) const {
    return data != x.data || pos != x.pos;
  }
  bool operator<(const block_map_iterator& x) const {
    return data == x.data && pos < x.pos;
  }
  

  };



template <typename Key, typename T, size_t N, 
  typename KeyMapFun, typename KeyUnmapFun>
class block_map {
protected:
  T data[N];
  const T empty_element;
  KeyMapFun key_map;
  KeyUnmapFun key_unmap;

public:
  
  typedef T mapped_type;
  typedef Key key_type;
  typedef std::pair<const key_type, mapped_type> value_type;
  
  typedef mapped_type* pointer;
  typedef const mapped_type* const_pointer;
  typedef mapped_type& reference;
  typedef const mapped_type& const_reference;
    
  typedef ptrdiff_t difference_type;
  typedef size_t size_type;
    
  typedef block_map_iterator<Key, T, KeyUnmapFun> iterator;
  typedef block_map_iterator<Key, T, KeyUnmapFun> const_iterator;
  
  typedef block_map<Key, T, N, KeyMapFun, KeyUnmapFun> base;
  

  
  iterator begin() { return iterator(data, 0, N, empty_element); }
    const_iterator begin() const {
      return const_iterator(data, 0, N, empty_element);
    }
  
    iterator end() { return iterator(data, N, N, empty_element); }
    const_iterator end() const {
      return const_iterator(data, N, N, empty_element);
    }
  

  
  size_type size() const {
  size_type s;
  for (size_type i = s = 0; i < N; ++i) 
    if (data[i] != empty_element) 
      ++s;
  return s;
  }
  
  size_type max_size() const { return N; }
  bool empty() const { return size() == 0; }
  

  
  reference operator[](key_type n) { return data[key_map(n)]; }
  const_reference operator[](key_type n) const {
    return data[key_map(n)];
  }
  
  iterator find(key_type n) {
    if (data[key_map(n)] != empty_element) {
      return iterator(data, key_map(n), N, empty_element);
    }
    return end();
  }
  
  size_type count(key_type n) { 
    return find(n) == end() ? 0 : 1; 
  }
  
  std::pair<iterator, iterator> equal_range(key_type n) {
    std::pair<iterator, iterator> ret;
    ret.first = ret.second = find(n);
    return ret;
  }
  

  
  bool operator==(const base& x) {
    return equal(data, &data[N], x.data);
  }
  
  bool operator!=(const base& x) { return !(*this == x); }
  
  bool operator<(base& x) {
    size_type s1 = size();
    size_type s2 = x.size();
    
    if (s1 < s2) return true;
    else if (s1 > s2) return false;
    else {
      iterator i = begin();
      iterator j = x.begin();
      for (; i < end(); ++i, ++j) {
        if (i->second < j->second) return true;
        else if (i->second > j->second) return false;
        }
      }
    return false;
  }
  
  bool operator>(base& x) { return x < (*this); }
  bool operator<=(base& x) { return !(x < (*this)); }
  bool operator>=(base& x) { return !((*this) < x); }
  

  
  std::pair<iterator, bool> insert(value_type p) {
  iterator iter(data, 0, N, empty_element);
  std::pair<iterator, bool> ret(iter, false);
  if (data[key_map(p.first)] == empty_element) {
    data[key_map(p.first)] = p.second;
    ret.first.pos = key_map(p.first);
    ret.second = true;
    }
  return ret;
  }
  
  template<typename InputIter>
  void insert(InputIter p, InputIter q) {
    while (p < q) insert(*p++);
  }
  
  void erase(key_type n) { data[key_map(n)] = empty_element; }
  void erase(iterator i) { data[key_map(i->first)] = empty_element; }
  void erase(iterator p, iterator q) {
    while (p < q) { data[key_map(p->first)] = empty_element; ++p; }
  }
  
  void swap(base& x) {
    swap_ranges(data, &data[N], x.data);
  }
  

  
  block_map(T e = T()) : empty_element(e) {
    for (size_type i = 0; i < N; ++i) data[i] = empty_element;
  }
  
  block_map(const base& x) : empty_element(x.empty_element) 
  { copy(x.data, &x.data[N], data); }
  
  template<typename InputIter>
  block_map(InputIter p, InputIter q, T e = T()) : 
  empty_element(e) {
    for (size_type i = 0; i < N; ++i) data[i] = empty_element;
    while (p < q) insert(*p++);
  }
  

};


#endif
