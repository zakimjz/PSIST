#ifndef __STREE_ARRAY__
#define __STREE_ARRAY__

#include "stree_config.h"
#include "block_map.h"
#include "stree.h"

struct dna_mapper : public std::unary_function<char, char> {
  dna_mapper() {}
  char operator()(char c) {
    switch(c) {
    case 'a':
    case 'A':
      return 0;
    case 'c':
    case 'C':
      return 1;
    case 'g':
    case 'G':
      return 2;
    case 't':
    case 'T':
      return 3;
    default:
      return 4;
    }
  }
};

struct dna_unmapper : public std::unary_function<char, char> {
  dna_unmapper() {}
  char operator()(char c) {
    switch(c) {
    case 0: return 'a';
    case 1: return 'c';
    case 2: return 'g';
    case 3: return 't';
    default: return '~';
    }
  }
};


template <typename Sequence, size_t AlphaSize,
  typename EqualKey, typename KeyMapFun, typename KeyUnmapFun>
struct _snode_array;

template <typename Sequence, size_t AlphaSize,
  typename EqualKey, typename KeyMapFun, typename KeyUnmapFun>
struct _sedge_array {
  typedef _snode_array<Sequence, AlphaSize, EqualKey, 
    KeyMapFun, KeyUnmapFun> node_type;
    
  typedef typename Sequence::difference_type difference_type;
  typedef long SequenceID;
  
  node_type* parent;
  node_type* child;
  SequenceID sequence_id;
  difference_type first, last;
  

  _sedge_array( difference_type in_first, difference_type in_last )
    : first(in_first), last(in_last) { }
};

template <typename Sequence, size_t AlphaSize,
  typename EqualKey, typename KeyMapFun, typename KeyUnmapFun>
struct _snode_array {
  typedef _snode_array<Sequence, AlphaSize, EqualKey, 
    KeyMapFun, KeyUnmapFun> node_type;
  
  typedef block_map<typename Sequence::value_type, 
    _sedge_array<Sequence, AlphaSize, EqualKey, KeyMapFun, KeyUnmapFun>*, 
    AlphaSize+1, KeyMapFun, KeyUnmapFun>
    edge_container_type;
  
  typedef typename Sequence::difference_type difference_type;
  typedef long SequenceID;
  
  edge_container_type edges;
  std::map<SequenceID, difference_type> labels;
  node_type* suffix_link;
  bool leaf;
  
  _snode_array() : suffix_link(NULL), leaf(false) { }
};

template <typename Sequence, 
  typename Sequence::value_type EndMarker, 
  size_t AlphaSize, typename KeyMapFun, typename KeyUnmapFun,
  typename EqualKey = std::equal_to<STREE_TYPENAME Sequence::value_type> >
class stree_array : 
  public _stree<Sequence, EndMarker, EqualKey,
    _snode_array<Sequence, AlphaSize, EqualKey, KeyMapFun, KeyUnmapFun>, 
    _sedge_array<Sequence, AlphaSize, EqualKey, KeyMapFun, KeyUnmapFun> >
{
private:
  typedef long SequenceID;
  typedef _snode_array<Sequence, AlphaSize, EqualKey, 
    KeyMapFun, KeyUnmapFun> node_type;
  typedef _sedge_array<Sequence, AlphaSize, EqualKey, 
    KeyMapFun, KeyUnmapFun> edge_type;
  typedef _stree<Sequence, EndMarker, EqualKey, 
    node_type, edge_type> base;

public:
  stree_array() : base() {}
  stree_array(const Sequence& s, SequenceID id) : base(s, id) {}
  ~stree_array() { clear(); }
};

#endif
