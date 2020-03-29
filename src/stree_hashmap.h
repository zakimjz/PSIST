#ifndef __STREE_HASH__
#define __STREE_HASH__

#include "stree_config.h"
#include <ext/hash_map>
#include "stree.h"

//#define HASHNS __gnu_cxx
using namespace __gnu_cxx;


template <typename Sequence, typename HashFun, typename EqualKey>
struct _snode_hashmap;

template <typename Sequence, typename HashFun, typename EqualKey>
struct _sedge_hashmap { 
  typedef _snode_hashmap<Sequence, HashFun, EqualKey> node_type;
  typedef typename Sequence::difference_type difference_type;
  typedef long SequenceID;
  
  node_type* parent;
  node_type* child;
  SequenceID sequence_id;
  difference_type first, last;
  
  _sedge_hashmap( difference_type in_first, difference_type in_last )
    : first(in_first), last(in_last) { }
};

template <typename Sequence, typename HashFun, typename EqualKey>
struct _snode_hashmap {
  typedef _snode_hashmap<Sequence, HashFun, EqualKey> node_type;
  typedef hash_map<typename Sequence::value_type, 
    _sedge_hashmap<Sequence, HashFun, EqualKey>*, HashFun, EqualKey>
    edge_container_type;
  
  typedef typename Sequence::difference_type difference_type;
  typedef long SequenceID;
  
  edge_container_type edges;
  std::map<SequenceID, difference_type> labels;
  node_type* suffix_link;
  bool leaf;
  
  _snode_hashmap() : suffix_link(NULL), leaf(false) { }
};

template <typename Sequence, 
  typename Sequence::value_type EndMarker, 
  typename EqualKey = std::equal_to<STREE_TYPENAME Sequence::value_type>,
  typename HashFun = hash<STREE_TYPENAME Sequence::value_type> >
class stree_hashmap : 
  public _stree<Sequence, EndMarker, EqualKey,
    _snode_hashmap<Sequence, HashFun, EqualKey>, 
    _sedge_hashmap<Sequence, HashFun, EqualKey> >
{
  private:
  typedef _snode_hashmap<Sequence, HashFun, EqualKey> node_type;
  typedef _sedge_hashmap<Sequence, HashFun, EqualKey> edge_type;
  typedef _stree<Sequence, EndMarker, EqualKey, 
  node_type, edge_type > base;
  
  public:
  typedef long SequenceID;
  stree_hashmap() : base() {}
  stree_hashmap(const Sequence& s, SequenceID id) : base(s, id) {}
  ~stree_hashmap() { clear(); }

};
#endif
