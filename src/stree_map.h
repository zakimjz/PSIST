#ifndef __STREE_MAP__
#define __STREE_MAP__

#include "stree_config.h"
#include <map>
#include "stree.h"


using namespace std;

template <typename Sequence, typename Compare, typename EqualKey>
struct _snode_map;

template <typename Sequence, typename Compare, typename EqualKey>
struct _sedge_map {
   typedef _snode_map<Sequence, Compare, EqualKey> node_type;
   typedef typename Sequence::difference_type difference_type;
   typedef long SequenceID;
  
   node_type* parent;
   node_type* child;
   SequenceID sequence_id;
   difference_type first, last;
  
   _sedge_map( difference_type in_first, difference_type in_last )
      : first(in_first), last(in_last) { }
   
   //added MJZ -- 9/18/04
  //friend ostream& operator << (ostream& fout, _sedge_map& edge);
};


//added MJZ -- 9/18/04
/*
template <typename Sequence, typename Compare, typename EqualKey>
ostream& operator << (ostream& fout,  
                      _sedge_map<Sequence, Compare, EqualKey>& edge)
{
   fout << edge.sequence_id << ", [" << edge.first << " " << edge.last << "] ";
   return fout;
}
*/


template <typename Sequence, typename Compare, typename EqualKey>
struct _snode_map {
   typedef _snode_map<Sequence, Compare, EqualKey> node_type;

   typedef std::map<typename Sequence::value_type, _sedge_map<Sequence,
                                                              Compare, EqualKey>*, Compare>
     edge_container_type; 
  
   typedef typename Sequence::difference_type difference_type;
   typedef long SequenceID;
  
   edge_container_type edges;
   std::map<SequenceID, difference_type> labels;
   std::vector<std::pair<difference_type,
                         SequenceID> > ilabels; //internal node
                                                //labels,including
                                                //labels of all children
   
   node_type* suffix_link;
   bool leaf;
  
   _snode_map() : suffix_link(NULL), leaf(false) { }
};


// //added MJZ -- 9/18/04
// template <typename Sequence, typename Compare, typename EqualKey>
// ostream& operator << (ostream& fout,  
//                       _snode_map<Sequence, Compare, EqualKey>& node)
// {
//    typedef _snode_map<Sequence, Compare, EqualKey> node_type;
//    typedef typename node_type::edge_container_type::iterator edge_iterator;
//    for (edge_iterator i = node.edges.begin(); i!=node.edges.end(); ++i) {
//       _sedge_map<Sequence, Compare, EqualKey>* eptr = i->second;
//       fout << i->first << " -- " << *eptr << endl;
//       cout << "\tCHILD " << *eptr->child << endl;
//    }
//    return fout;
// }



template <typename Sequence, 
          typename EndMarker, 
	  typename EqualKey = std::equal_to<STREE_TYPENAME 
Sequence::value_type>,
          typename Compare = std::less<STREE_TYPENAME Sequence::value_type> >
class stree_map : 
   public _stree<Sequence, EndMarker, EqualKey,
                 _snode_map<Sequence, Compare, EqualKey>, 
                 _sedge_map<Sequence, Compare, EqualKey> >
{
  private:
  
   typedef _snode_map<Sequence, Compare, EqualKey> node_type;
   typedef _sedge_map<Sequence, Compare, EqualKey> edge_type;
   typedef _stree<Sequence, EndMarker, EqualKey, 
                  node_type, edge_type > base;



public:
  typedef long SequenceID;
  stree_map() : base() {}
  stree_map(EndMarker em) : base(em) {}
  stree_map(const Sequence& s, SequenceID id) : base(s, id) {}
  ~stree_map() { 
    //clear();
  }
};


#endif
