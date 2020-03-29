
#ifndef __STREE__
#define __STREE__

#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include "stree_structures.h"

using namespace std;

template <typename Sequence, 
          typename EndMarker, 
          typename EqualKey,
          typename Node,
          typename Edge> 
 
class _stree { 
 protected: 
  
   typedef typename Sequence::value_type value_type;
   typedef typename Sequence::difference_type difference_type;
   typedef typename Sequence::size_type size_type;
   typedef long SequenceID;
   typedef typename Node::edge_container_type::iterator edge_iterator;


   //--------------moved to public-------------------
   //Node* root;
   //std::map<SequenceID, Sequence> sequences;
   //std::map<SequenceID, size_type> lengths;
   //--------------moved to public-------------------
   
  
   void linear_time_construction(SequenceID id); 
   void split_edge( Edge&, difference_type );
   Node* find_node( const Sequence&);
  

  
   difference_type edge_length( Edge& edge ) {
      if (edge.child->leaf) {
         return lengths[edge.sequence_id]-edge.first+1;
      } else {
         return edge.last-edge.first+1;
      }
   }
  
 
  
   void add_leaf( Node& node, difference_type index, 
                  difference_type label, SequenceID id ) {
      Edge* eptr = node.edges[(sequences[id])[index]] = 
         new Edge(index, index);
      eptr->sequence_id = id;
      eptr->parent = &node;
      eptr->child = new Node();
      eptr->child->labels[id] = label;
      eptr->child->leaf = true;
   }
  
 
  
   void find_leaves( Node& node, std::vector<std::pair<difference_type,
                     SequenceID> >& matches)
   {
      if (node.leaf) {
         for (typename std::map<SequenceID, difference_type>::
                 iterator i = node.labels.begin(); i != node.labels.end(); ++i) 
            {
               matches.push_back(std::make_pair(i->second, i->first));
            }
      } else {
         for(edge_iterator i = node.edges.begin(); i != node.edges.end(); ++i) 
            {
               find_leaves(*(i->second->child), matches);
            }
      }
   }
  

  
   void find_leaves_in( Node& node, std::vector<difference_type>& matches, 
                        SequenceID id)
   {
      if (node.leaf) {
         typename std::map<SequenceID, difference_type>::iterator i = node.labels.find(id);
         if (i != node.labels.end()) matches.push_back(i->second);
      } else {
         for(edge_iterator i = node.edges.begin(); i != node.edges.end(); ++i) 
            {
               find_leaves_in(*(i->second->child), matches, id);
            }
      }
   }
  
 
  
   void init(SequenceID id) {
      root = new Node();
      root->suffix_link = root;
      root->labels[id] = 0;
   }
  
   
   void clear(Node* nptr);
   void fix_label(Edge*, SequenceID);
   void erase(Node* nptr, SequenceID id);
   bool erase_leaf(Edge* eptr, SequenceID id);
   void join_edge(Edge* eptr, SequenceID id);
   
   
 public:

   
   Node* root;
   std::map<SequenceID, Sequence> sequences;
   std::map<SequenceID, size_type> lengths;
   
   EndMarker  emarker;
   
   _stree() : root(NULL){}
 
   _stree(EndMarker& em) : root(NULL),emarker(em) { }
   
   _stree(const 
          _stree<Sequence, EndMarker,  EqualKey, Node, Edge>
          & s) {
      if (&s != this) {
         root = NULL;
         for (typename std::map<SequenceID, Sequence>::const_iterator 
                 i = s.sequences.begin();
              i != s.sequences.end(); ++i)
	   insert(i->second, i->first);
	 emarker = s.emarker;
      }
   }  
   
   _stree(const Sequence& s, SequenceID id) {
      root = NULL;
      insert(s, id);
   }


   void insert(const Sequence& s, SequenceID id) {
      sequences[id] = s;
      if (*(s.end()-1) != emarker)
         sequences[id].insert(sequences[id].end(),emarker );
      if (root == NULL) init(id);
      linear_time_construction(id);
   }
  

  
   void clear() {
      if (root == NULL) return;
      clear(root);
      sequences.clear();
      lengths.clear();
      root = NULL;
   }
  

  
   void erase(SequenceID id) {
      if (sequences.find(id) == sequences.end()) 
         return;
      if (sequences.size() == 1) 
         clear();
      else {
         erase(root, id);
         sequences.erase(id);
         lengths.erase(id);
      }
   }
  

  
   bool empty() const { return root == NULL; }
   size_type size() const { 
      size_type s = 0;
      for (typename std::map<SequenceID, size_type>::const_iterator i = lengths.begin();
           i != lengths.end(); ++i)
         {
            s += i->second;
         }
      return s;
   }
   size_type max_size() const {
      Sequence s;
      return s.max_size();
   }
   Node* get_root() { return root; }
   const Sequence& get_sequence(SequenceID n) { return sequences[n]; }
  

  
   std::vector< std::pair<difference_type, SequenceID> > 
   find( const Sequence& pattern) {
      std::vector< std::pair<difference_type, SequenceID> > matches;
      Node* nptr = find_node(pattern);
      if (nptr != NULL) find_leaves(*nptr, matches);
      return matches;
   }
  

   // exact matching

   void find_lcs( const Sequence& pattern, size_type lcs_min, std::vector< _matches >& matches, value_type pref) {
     
       //std::vector< std::pair<difference_type, SequenceID> > matches;
       //vector<_matches > matches;
    
       //Node* nptr = find_node(pattern);
       const size_type p_size(pattern.size());
       Node* nptr = root;
       size_type i(0);
       size_type elength(0);
       edge_iterator iter;
       EqualKey equal;
       
       while ( i < lcs_min &&
	      (iter = nptr->edges.find(pattern[i])) != nptr->edges.end()) {
	   Edge* eptr = iter->second;
	   size_type e_max = eptr->first + edge_length(*eptr); 
	   size_type m = eptr->first + 1;
	   ++i;
	   
	   //elength = (e_max - m) - (p_size -i);
	   if ( m < e_max) {
	       const Sequence& s(sequences[eptr->sequence_id]);
	       while ( m < e_max && equal(s[m], pattern[i])) {
		   ++i; ++m;
	       }
	   }
	   
	   nptr = eptr->child; 
       }

       if(nptr == NULL || i< lcs_min)  return;  //empty
       //cout << "--i--" << i << endl;

       while( i < p_size && nptr != NULL){
	   Node* grownode;
	   size_type elength = i;

	   //cout << "--i--" << i << endl;
	   //bool morematch = false;
	  
	   if((iter = nptr->edges.find(pattern[i])) != nptr->edges.end()){
	       
	     Edge* eptr = iter->second;
	     size_type e_max = eptr->first + edge_length(*eptr); 
	     size_type m = eptr->first + 1;
	     ++i;
	     
	     if ( m < e_max) {
	       const Sequence& s(sequences[eptr->sequence_id]);
	       while ( m < e_max && equal(s[m], pattern[i])) {
		 ++i; ++m;
	       }
	     }
	     
	     grownode = eptr->child;
	     
	     //find the leaves for all the other unmatching nodes
	     for(edge_iterator iter2 = nptr->edges.begin(); iter2 != nptr->edges.end(); ++iter2)
	       if(iter2 != iter )
		 find_lcs_leaves( *(iter2->second->child), matches, elength, pref);
	     
	   }
	   else{ //no more matches
	     if(nptr->leaf)
	       find_lcs_leaves( *nptr, matches, i, pref);
	     break;
	   }
	   
	   //find_lcs_leaves(*(iter->second->child), matches, elength);
	   //}

	   nptr = grownode;	  
	   
	   //if(morematch == false){
	   //  if(nptr->leaf)
	   //   find_lcs_leaves(*nptr, matches, i);
	   //  break;
	   //}
       
       }
       
       //cout << "--i--" << i << endl;
       if(i == p_size)
	   find_lcs_leaves(*nptr, matches, i, pref);
       
       //return matches;
   }
   



   void find_lcs_leaves(Node& node, std::vector<_matches >& matches, size_type elength, value_type pref)
  {
     // cout << "m size --" << matches.size() << "\t" << elength <<endl;
     if (node.leaf) {
       for (typename std::map<SequenceID, difference_type>::
	      iterator it = node.labels.begin(); it!= node.labels.end(); ++it) 
	 {
	   
	   long ref_prefix = it->second - 1;
	   if( ref_prefix >= 0){
	     Sequence& s(sequences[it->first]);
	     EqualKey equal;
	     if (equal(pref, s[ref_prefix]) )
	       continue;
	   }
	   
	   _matches m;
	   m.lable  = it->second;
	   m.length = elength;
	   m.sequenceid = it->first;
	   //cout << "--add---" << m.lable << " " << m.length << " " << m.sequenceid << endl;
	   matches.push_back(m);
	 }
     } else {
       for(edge_iterator it = node.edges.begin(); it!= node.edges.end(); ++it) 
	 {
	   //elength = edge_length(i->second);
	   find_lcs_leaves(*(it->second->child), matches, elength, pref);
	 }
     }
   }
   
   
   
   void find_appr_lcs_leaves(Node& node, std::vector<_matches >& matches, size_type elength, float epsilon, value_type pref)
  {
    // cout << "m size --" << matches.size() << "\t" << elength <<endl;
    if (node.leaf) {
       for (typename std::map<SequenceID, difference_type>::
	      iterator it = node.labels.begin(); it!= node.labels.end(); ++it) 
	 {
	   
	   long ref_prefix = it->second - 1;
	   if( ref_prefix >= 0){
	     Sequence& s(sequences[it->first]);
	     EqualKey equal(epsilon);
	     if (equal(pref, s[ref_prefix]) )
	       continue;
	   }
	   
	   _matches m;
	   m.lable  = it->second;
	   m.length = elength;
	   m.sequenceid = it->first;
	   //cout << "--add---" << m.lable << " " << m.length << " " << m.sequenceid << endl;
	   matches.push_back(m);
	 }
     } else {
       for(edge_iterator it = node.edges.begin(); it!= node.edges.end(); ++it) 
	 {
	   //elength = edge_length(i->second);
	   find_appr_lcs_leaves(*(it->second->child), matches, elength, epsilon, pref);
	 }
     }
   }
   
   
   void find_appr_lcs_node(const Sequence& pattern, size_type lcs_min, Node* nptr, std::vector<_matches>& matches, size_type i, float& epsilon, value_type pref) {
     const size_type p_size(pattern.size());
     
     
     if(i >= p_size){
       find_appr_lcs_leaves(*nptr, matches, i, epsilon, pref);
       return;
     }
     
     if (nptr->leaf) {
       if( i >= lcs_min )
	 find_appr_lcs_leaves(*nptr, matches, i, epsilon, pref);
     }
     else { // its a node
       
       EqualKey equal(epsilon); //same to all function, as a parameter.
       
       for(edge_iterator iter = nptr->edges.begin(); iter != nptr->edges.end(); iter++){
	
	 if( equal(pattern[i], iter->first) == true ){
	   size_type elength = i;
	   Edge* eptr = iter->second;
	   size_type e_max = eptr->first + edge_length(*eptr); 
	   size_type m = eptr->first + 1;
	   ++elength;
	   
	   if ( m < e_max) {
	     const Sequence& s(sequences[eptr->sequence_id]);
	     while ( m < e_max && equal(s[m], pattern[elength])) {
	       ++elength; ++m;
	     }
	   }
	   
	   //nptr = eptr->child;
	   find_appr_lcs_node(pattern, lcs_min, eptr->child, matches, elength, epsilon, pref);
	 }
	 else{
	   if( i >= lcs_min )
	     find_appr_lcs_leaves(*nptr, matches, i, epsilon, pref); 
	 }
       }// end for
       
     }// end node
   }
   
   
   
   // inexact matching
   
   void find_appr_lcs( const Sequence& pattern, size_type lcs_min,  std::vector< _matches >& matches, float epsilon, value_type pref) {
     
     find_appr_lcs_node(pattern, lcs_min, root, matches, 0, epsilon, pref);
   }
   
   
   void find_appr_lcs_leaves_stree( Node& node, map<SequenceID, Sequence>& qseq, vector<pair<difference_type, SequenceID> >& qlables, std::vector<_matches >& matches, size_type elength, float epsilon =0.0)
   {
     // cout << "m size --" << matches.size() << "\t" << elength <<endl;
     if (node.leaf) {
       for (typename std::map<SequenceID, difference_type>::
	      iterator it = node.labels.begin(); it!= node.labels.end(); ++it) 
	 {
	   
	   int ql_size = qlables.size();
	   for( int ql=0; ql < ql_size; ++ql){
	     
	     SequenceID ref_prefix = it->second - 1;
	     SequenceID q_prefix = qlables[ql].first - 1;

	     if( ref_prefix >= 0 &&  q_prefix >=0 ){
	       Sequence& s(sequences[it->first]);
	       Sequence& qs(qseq[qlables[ql].second]);
	       
	       EqualKey equal(epsilon);
	       if (equal(qs[q_prefix], s[ref_prefix]) )
		 continue;
	     }
	     
	     _matches m;
	     m.lable  = it->second;
	     m.length = elength;
	     m.sequenceid = it->first;
	     m.querylable = qlables[ql].first;
	     m.queryid = qlables[ql].second;
	     matches.push_back(m);
	   }
	 }
     }
     else
       for(edge_iterator it = node.edges.begin(); it!= node.edges.end(); ++it) 
	 find_appr_lcs_leaves_stree(*(it->second->child), qseq, qlables, matches, elength, epsilon);
   }
   
   
   void find_appr_lcs_stree(Node* qnode, map<SequenceID, Sequence>& qseq, size_type lcs_min,  std::vector<_matches >& matches, float epsilon) {

     for(edge_iterator qiter = qnode->edges.begin(); qiter != qnode->edges.end(); qiter++){
       for(edge_iterator iter = root->edges.begin(); iter != root->edges.end(); iter++){
	 Edge* qeptr = qiter->second;
	 Edge* eptr  = iter->second;
	 find_appr_lcs_node_stree(qeptr,qeptr->first,eptr,eptr->first,0, qseq, lcs_min, matches,epsilon);
       }
     }
     
   }
   
   
   
   void find_appr_lcs_node_stree(Edge* qedge, size_type qi, Edge* edge, size_type pi, size_type i, map<SequenceID, Sequence>& qseq, size_type lcs_min, std::vector<_matches>& matches, float& epsilon) {
     
     Node* qnode = qedge->child;
     Node* node = edge->child;
     EqualKey equal(epsilon); //same to all function, as a parameter.
     
     const Sequence& qs(qseq[qedge->sequence_id]);
     const Sequence& s(sequences[edge->sequence_id]);

     if( (equal(qs[qi], s[pi]) == true) ) {
       
       //cout << "qs" << qs[qi] << endl;
       size_type qe_max = qedge->first + edge_length(*qedge); 
       size_type e_max = edge->first + edge_length(*edge); 
       
       size_type qm = qi + 1;
       size_type m  = pi + 1;
       ++i;
       
       if ( m < e_max && qm < qe_max) {
	 while ( m < e_max && qm < qe_max && equal(qs[qm], s[m])) {
	   ++qm; ++m; ++i;
	 }
       }
       
       //i += (qm - qi);
       
       if( m == e_max && qm == qe_max ){

	 //if(qnode->leaf && node->leaf){ //both leaves, (EndMarker) -1 == -1, 
	 //if(i >= lcs_min){            //which makes length (i) is +1
	 //  vector<pair<difference_type, SequenceID> > qlables;
	 //  find_leaves(*qnode, qlables);
	 //find_appr_lcs_leaves_stree(*node, qseq,qlables,matches,i-1,epsilon);
	 //}
	 //return;
	 //}
	 
	 if(qnode->leaf || node->leaf){
	   if(i >= lcs_min){
	     vector<pair<difference_type, SequenceID> > qlables;
	     find_leaves(*qnode, qlables);
	     find_appr_lcs_leaves_stree(*node,qseq, qlables,matches,i-1,epsilon);
	   }
	 }
	 else
	   for(edge_iterator qiter = qnode->edges.begin(); qiter != qnode->edges.end(); qiter++)
	     for(edge_iterator iter = node->edges.begin(); iter != node->edges.end(); iter++)
	       find_appr_lcs_node_stree(qiter->second, qiter->second->first,iter->second,iter->second->first, i, qseq, lcs_min, matches,epsilon);
	 
	 return;
       }
       
       if( m == e_max ){
	 
	 if(node->leaf){
	   if(i >= lcs_min){
	     vector<pair<difference_type, SequenceID> > qlables;
	     find_leaves(*qnode, qlables);
	     find_appr_lcs_leaves_stree(*node,qseq, qlables,matches,i-1,epsilon);
	   }
	 }
	 else
	   for(edge_iterator iter = node->edges.begin(); iter != node->edges.end(); iter++)
	     find_appr_lcs_node_stree(qedge, qm,iter->second,iter->second->first, i, qseq, lcs_min, matches,epsilon);
	 return;
       }
       
       
       if(qm == qe_max ){
	 
	 if(qnode->leaf){
	   if(i >= lcs_min){
	     vector<pair<difference_type, SequenceID> > qlables;
	     find_leaves(*qnode, qlables);
	     find_appr_lcs_leaves_stree(*node,qseq, qlables,matches,i-1,epsilon);
	   }
	 }
	 else
	   for(edge_iterator qiter = qnode->edges.begin(); qiter != qnode->edges.end(); qiter++)
	     find_appr_lcs_node_stree(qiter->second, qiter->second->first,edge,m, i, qseq, lcs_min, matches,epsilon);
	 return;
       }
       
       
       if( m != e_max && qm != qe_max ){
	 
	 if(i >= lcs_min){
	   vector<pair<difference_type, SequenceID> > qlables;
	   find_leaves(*(qedge->child), qlables);
	   find_appr_lcs_leaves_stree(*(edge->child),qseq, qlables,matches,i,epsilon);
	 }
	 return;
       }
       
     }
     else{
       if(i >= lcs_min){
	 vector<pair<difference_type, SequenceID> > qlables;
	 find_leaves(*qnode, qlables);
	 find_appr_lcs_leaves_stree(*node,qseq, qlables,matches,i,epsilon);
       }       
     }
     

   }
   
   
   
   
   
   void find_lcs_leaves_stree( Node& node, map<SequenceID, Sequence>& qseq, vector<pair<difference_type, SequenceID> >& qlables, std::vector<_matches >& matches, size_type elength)
   {
     // cout << "m size --" << matches.size() << "\t" << elength <<endl;
     if (node.leaf) {
       for (typename std::map<SequenceID, difference_type>::
	      iterator it = node.labels.begin(); it!= node.labels.end(); ++it) 
	 {
	   
	   int ql_size = qlables.size();
	   for( int ql=0; ql < ql_size; ++ql){
	     
	     SequenceID ref_prefix = it->second - 1;
	     SequenceID q_prefix = qlables[ql].first - 1;

	     if( ref_prefix >= 0 &&  q_prefix >=0 ){
	       Sequence& s(sequences[it->first]);
	       Sequence& qs(qseq[qlables[ql].second]);
	       
	       EqualKey equal;
	       if (equal(qs[q_prefix], s[ref_prefix]) )
		 continue;
	     }
	     
	     _matches m;
	     m.lable  = it->second;
	     m.length = elength;
	     m.sequenceid = it->first;
	     m.querylable = qlables[ql].first;
	     m.queryid = qlables[ql].second;
	     matches.push_back(m);
	   }
	 }
     }
     else
       for(edge_iterator it = node.edges.begin(); it!= node.edges.end(); ++it) 
	 find_lcs_leaves_stree(*(it->second->child), qseq, qlables, matches, elength);
   }
   
   
   
   void find_lcs_stree(Node* qnode, map<SequenceID, Sequence>& qseq, size_type lcs_min,  std::vector<_matches >& matches) {

     for(edge_iterator qiter = qnode->edges.begin(); qiter != qnode->edges.end(); qiter++){
       Edge* qeptr = qiter->second;
       find_lcs_node_stree(qeptr,qeptr->first, root, 0, 0, qseq, lcs_min, matches);
     }
     
   }
   
   
   
   void find_lcs_node_stree(Edge* qedge, size_type qi, Node* nptr, size_type pi, size_type i, map<SequenceID, Sequence>& qseq, size_type lcs_min, std::vector<_matches>& matches) {
     
          
     EqualKey equal; //same to all function, as a parameter.
     edge_iterator iter;
     
     const Sequence& qs(qseq[qedge->sequence_id]);
     //const Sequence& s(sequences[edge->sequence_id]);
     
     
     if( (iter = nptr->edges.find(qs[qi])) != nptr->edges.end() ) {
       
       //cout << "qs" << qs[qi] << endl;

       int elength = i;
       
       Node* qnode = qedge->child;
       size_type qe_max = qedge->first + edge_length(*qedge); 
       
       Edge* eptr = iter->second;
       Node* node = eptr->child;
       size_type e_max = eptr->first + edge_length(*eptr); 
       //size_type m = pi + 1;
       size_type m = eptr->first + 1;
       
       //size_type e_max = edge->first + edge_length(*edge); 
       
       size_type qm = qi + 1;
       //size_type m  = pi + 1;
       ++i;
       
       const Sequence& s(sequences[eptr->sequence_id]);
       
       if ( m < e_max && qm < qe_max) {
	 while ( m < e_max && qm < qe_max && equal(qs[qm], s[m])) {
	   ++qm; ++m; ++i;
	 }
       }
       
       //i += (qm - qi);
       
       if( m == e_max && qm == qe_max ){
	 
	 //if(qnode->leaf && node->leaf){ //both leaves, (EndMarker) -1 == -1, 
	 //if(i >= lcs_min){            //which makes length (i) is +1
	 //  vector<pair<difference_type, SequenceID> > qlables;
	 //  find_leaves(*qnode, qlables);
	 //  find_lcs_leaves_stree(*node, qseq,qlables,matches,i-1,0);
	 //}
	 //return;
	 //}
	 
	 if(qnode->leaf || node->leaf){
	   if(i >= lcs_min){
	     vector<pair<difference_type, SequenceID> > qlables;
	     find_leaves(*qnode, qlables);
	     find_lcs_leaves_stree(*node,qseq, qlables,matches,i-1);
	   }
	 }
	 else{

	   for(edge_iterator qiter = qnode->edges.begin(); qiter != qnode->edges.end(); qiter++)
	     find_lcs_node_stree(qiter->second, qiter->second->first, node, m, i, qseq, lcs_min, matches);
	 }
	 
	 return;
       }
       
       if( m == e_max ){
	 
	 if(node->leaf){
	   if(i >= lcs_min){
	     vector<pair<difference_type, SequenceID> > qlables;
	     find_leaves(*qnode, qlables);
	     find_lcs_leaves_stree(*node,qseq, qlables,matches,i-1);
	   }
	 }
	 else
	   find_lcs_node_stree(qedge, qm, node,m, i, qseq, lcs_min, matches);
	 return;
       }
       
       
       if(qm == qe_max ){
	 
	 if(qnode->leaf){
	   if(i >= lcs_min){
	     vector<pair<difference_type, SequenceID> > qlables;
	     find_leaves(*qnode, qlables);
	     find_lcs_leaves_stree(*node,qseq, qlables,matches,i-1);
	   }
	 }
	 else{
	   edge_iterator qiter;
	   while( (qiter = qnode->edges.find(s[m])) != qnode->edges.end()){

	     if(i >=lcs_min )
	       for(edge_iterator qiter2 = qnode->edges.begin(); qiter2 != qnode->edges.end(); ++qiter2)
		 if(qiter2 != qiter ){
		   
		   vector<pair<difference_type, SequenceID> > qlables;
		   find_leaves(*(qiter2->second->child), qlables);
		   find_lcs_leaves_stree( *node, qseq, qlables, matches, i);
		 }
       
	     
	     
	     Edge* qeptr = qiter->second;
	     size_type ne_max = qeptr->first + edge_length(*qeptr); 
	     size_type nm = qeptr->first + 1;
	     ++i;
	     ++m;
	     	     
	     //elength = (e_max - m) - (p_size -i);
	     if ( m < e_max && nm < ne_max) {
	       const Sequence& qs(qseq[qeptr->sequence_id]);
	       while ( m < e_max &&  nm < ne_max && equal(s[m], qs[nm])) {
		 ++i; ++nm; ++m;
	       }
	     }
	     
	     qnode = qeptr->child; 

	     if( m < e_max && nm == ne_max){
	       if(qnode->leaf) {
		 if(i >= lcs_min){
		   vector<pair<difference_type, SequenceID> > qlables;
		   find_leaves(*qnode, qlables);
		   find_lcs_leaves_stree(*node,qseq, qlables,matches,i-1);
		   return;
		 }
	       }
	       else  
		 continue;
	     }
	     
	     if( m == e_max && nm == ne_max ){
	       
	       if(qnode->leaf || node->leaf) {
		 if(i >= lcs_min){
		   vector<pair<difference_type, SequenceID> > qlables;
		   find_leaves(*qnode, qlables);
		   find_lcs_leaves_stree(*node,qseq, qlables,matches,i-1);
		 }
	       }
	       else
		 for(qiter = qnode->edges.begin(); qiter != qnode->edges.end(); qiter++)
		   find_lcs_node_stree(qiter->second, qiter->second->first, node,m,i, qseq, lcs_min, matches);
	       return;
	     }

	     if( m == e_max && nm < ne_max ){
	       
	       if(node->leaf){
		 if(i >= lcs_min){
		   vector<pair<difference_type, SequenceID> > qlables;
		   find_leaves(*qnode, qlables);
		   find_lcs_leaves_stree(*node,qseq, qlables,matches,i-1);
		 }
	       }
	       else
		 find_lcs_node_stree(qeptr, nm, node, m,i, qseq, lcs_min, matches);
	       return;
	     }


	     if( m < e_max && qm < qe_max ){
	       
	       if(i >= lcs_min){
		 vector<pair<difference_type, SequenceID> > qlables;
		 find_leaves(*qnode, qlables);
		 find_lcs_leaves_stree(*node,qseq, qlables,matches,i);
	       }
	       return;
	     }
	     

	   }// end while


	   if(i >= lcs_min){
	     vector<pair<difference_type, SequenceID> > qlables;
	     find_leaves(*qnode, qlables);
	     find_lcs_leaves_stree(*node,qseq, qlables,matches,i);
	   }
	 
	 }
	 
	 return;
       }
       
       
       if( m != e_max && qm != qe_max ){
	 
	 if(i >= lcs_min){
	   vector<pair<difference_type, SequenceID> > qlables;
	   find_leaves(*qnode, qlables);
	   find_lcs_leaves_stree(*node,qseq, qlables,matches,i);
	 }
	 return;
       }
       
       if(i >=lcs_min )
	 for(edge_iterator iter2 = nptr->edges.begin(); iter2 != nptr->edges.end(); ++iter2)
	   if(iter2 != iter ){
	     
	     vector<pair<difference_type, SequenceID> > qlables;
	     find_leaves(*qnode, qlables);
	     find_lcs_leaves_stree( *(iter2->second->child), qseq, qlables, matches, elength);
	   }
       
	     
     }
     else{
       if(i >= lcs_min){
	 vector<pair<difference_type, SequenceID> > qlables;
	 find_leaves(*(qedge->child), qlables);
	 find_lcs_leaves_stree(*nptr,qseq, qlables,matches,i);
       }       
     }
     

   }
   
   
   




   std::vector<difference_type> 
   find_in(const Sequence& pattern, SequenceID id) {
      std::vector<difference_type> matches;
      Node* nptr = find_node(pattern);
      if (nptr != NULL) find_leaves_in(*nptr, matches, id);
      return matches;
   }
  

  
   void swap( 
             _stree<Sequence, EndMarker,  EqualKey, Node, Edge>& s ) {
      Node* temp_root = s.root;
      std::map<SequenceID, Sequence> temp_sequences = s.sequences;
      std::map<SequenceID, size_type> temp_lengths = s.lengths;
      
      s.root = root;
      s.sequences = sequences;
      s.lengths = lengths;
      
      root = temp_root;
      sequences = temp_sequences;
      lengths = temp_lengths;
   }
  

   //added by MJZ -- 9/18/04
   void print(Node &nptr, int depth=0);
   
   friend ostream& operator << (ostream& os, 
                                _stree<Sequence, 
                                EndMarker,  EqualKey, Node, Edge>& t)
       {
	   t.print(*t.root, 0);
	   return os;
       }
   
   
   void propagate_labels(Node &nptr);
};


template <typename Sequence, 
          typename EndMarker, 

          typename EqualKey,
          typename Node,
          typename Edge>
void _stree<Sequence, EndMarker,  EqualKey, Node, Edge>
::print(Node &nptr, int depth)
{
   if (nptr.leaf) {
      for (typename std::map<SequenceID, difference_type>::
              iterator i = nptr.labels.begin(); i != nptr.labels.end(); ++i){
         for (int d=0; d < depth; ++d) cout << "\t";
         cout << "leaf: " << i->second << " " << i->first << endl;
      }
   } 
   else{
      for (int d=0; d < depth; ++d) cout << "\t";
      cout << "lbls: ";
      for (unsigned int j=0; j < nptr.ilabels.size(); ++j){
	  //cout<<"===" <<endl;
         cout << nptr.ilabels[j].first << " " 
              << nptr.ilabels[j].second << ", ";
      }
      cout << endl;
         
      for (edge_iterator i = nptr.edges.begin(); i!=nptr.edges.end(); ++i) {
         for (int d=0; d < depth; ++d) cout << "\t";
         cout << i->first << " -- " << *(i->second) << endl;
         print(*(i->second->child), depth+1);
      }
   }
}

template <typename Sequence, 
          typename EndMarker, 

          typename EqualKey,
          typename Node,
          typename Edge>
void _stree<Sequence, EndMarker,  EqualKey, Node, Edge>
::propagate_labels(Node &nptr)
{
   if (nptr.leaf) return; //we do nothing for leaf nodes

   //first find all child labels
   for (edge_iterator i = nptr.edges.begin(); i!=nptr.edges.end(); ++i) {
      if (!i->second->child->leaf) propagate_labels(*(i->second->child));
   }
   
   //now add all child labels to self label

   for (edge_iterator i = nptr.edges.begin(); i!=nptr.edges.end(); ++i) {
      Node *child = i->second->child;
      if (child->leaf){
         for (typename std::map<SequenceID, difference_type>::
                 iterator j = child->labels.begin(); 
              j != child->labels.end(); ++j){
            nptr.ilabels.push_back(std::make_pair(j->second, j->first));
         }
      }
      else{
         for (unsigned int j = 0; j < child->ilabels.size(); ++j){
            nptr.ilabels.push_back(child->ilabels[j]);
         }         
      }
   } 
}


template <typename Sequence, 
          typename EndMarker, 

          typename EqualKey,
          typename Node,
          typename Edge> 

Node* 
_stree<Sequence, EndMarker, EqualKey, Node, Edge>
::find_node( 
            const Sequence& pattern) {
   const size_type p_size(pattern.size());
   Node* nptr = root;
   size_type i(0);
   edge_iterator iter;
   EqualKey equal;
  
   while (i < p_size && 
          (iter = nptr->edges.find(pattern[i])) != nptr->edges.end()) {
      Edge* eptr = iter->second;
      size_type e_max = eptr->first + edge_length(*eptr); 
      size_type m = eptr->first + 1;
      ++i;

      if (i < p_size && m < e_max) {
         const Sequence& s(sequences[eptr->sequence_id]);
         while (i < p_size && m < e_max && equal(s[m], pattern[i])) {
            ++i; ++m;
         }
      }

      nptr = eptr->child; 
   }
   if (i == p_size) 
      return nptr;
   return NULL;
}


/*

template <typename Sequence, 
          typename EndMarker, 
          typename EqualKey,
          typename Node,
          typename Edge> 

Node* 
_stree<Sequence, EndMarker, EqualKey, Node, Edge>
::find_lcs_node( 
            const Sequence& pattern) {
   const size_type p_size(pattern.size());
   Node* nptr = root;
   size_type i(0);
   edge_iterator iter;
   EqualKey equal;
  
   while (i < p_size && 
          (iter = nptr->edges.find(pattern[i])) != nptr->edges.end()) {
      Edge* eptr = iter->second;
      size_type e_max = eptr->first + edge_length(*eptr); 
      size_type m = eptr->first + 1;
      ++i;

      if (i < p_size && m < e_max) {
         const Sequence& s(sequences[eptr->sequence_id]);
         while (i < p_size && m < e_max && equal(s[m], pattern[i])) {
            ++i; ++m;
         }
      }

      nptr = eptr->child; 
   }
   if (i == p_size) 
      return nptr;
   return NULL;
}

*/


template <typename Sequence, 
          typename EndMarker, 

          typename EqualKey,
          typename Node,
          typename Edge> 

void
_stree<Sequence, EndMarker,  EqualKey, Node, Edge>::
split_edge( Edge & oldedge,difference_type splitindex)
{
   difference_type index = splitindex + 1;
  
   Node * splitnode = new Node();
   Edge * splitedge;
   if (oldedge.child->leaf)  { 
      splitedge = new Edge(index, lengths[oldedge.sequence_id]);
   } else {
      splitedge = new Edge(index, oldedge.last);
   }
   splitedge->sequence_id = oldedge.sequence_id;
   splitnode->edges[(sequences[oldedge.sequence_id])[index]]=splitedge;
   splitedge->child = oldedge.child;
   splitedge->parent = splitnode;
  
   oldedge.child = splitnode;
   oldedge.last  = splitindex;
 
   splitnode->suffix_link = NULL;
}




template <typename Sequence, 
          typename EndMarker, 

          typename EqualKey,
          typename Node,
          typename Edge> 

void _stree<Sequence, EndMarker,  EqualKey, Node, Edge>
::linear_time_construction(SequenceID id) {
  
   size_type i;
   size_type j;
   size_type m;
   size_type beta;
   size_type jumps;
   bool SEARCHING, END_EXPANSION_PHASE;
  
   EqualKey equal;
   
   Sequence &cur_sequence = sequences[id];
   size_type length = cur_sequence.size()-1;
   
   size_type elen;
  
   Node* nptr;
   Edge* eptr;
   Node* last_suffix;
   
   for (i=0, j=0, nptr=root, last_suffix = NULL, jumps = 0; 
        i <=length; ++i) { 
       END_EXPANSION_PHASE = false;
       lengths[id] = i;
      
    
       while (!END_EXPANSION_PHASE) {
	  
	   m = j+jumps;
	   beta = i-m; 
	   SEARCHING = true;
	   while (beta != 0 && SEARCHING) {
	       eptr = nptr->edges[cur_sequence[m]];
	       elen = edge_length(*eptr);
	       if (beta >= elen) {
		   nptr = eptr->child;
		   beta = beta-elen;
		   m = m+elen;
		   jumps = jumps+elen;
	       } else {
		   SEARCHING = false;
	       }
	   } 
      
    
      
	   if (nptr->edges.find(cur_sequence[m]) != nptr->edges.end()) {
	       eptr = nptr->edges[cur_sequence[m]];
            if (equal(cur_sequence[i], 
                      (sequences[eptr->sequence_id])[eptr->first + beta])) {
          
               if (i != length) { END_EXPANSION_PHASE = true; }
               if (beta == 0) {
                  if (last_suffix != NULL) {
                     last_suffix->suffix_link = nptr;
                     last_suffix = NULL;
                  }
               }
               if (i == length) { eptr->child->labels[id] = j; }
          
      
            } else {
          
               split_edge(*eptr, eptr->first+beta-1); 
               add_leaf(*eptr->child, i, j, id);
               if (last_suffix != NULL) {
                  last_suffix->suffix_link = eptr->child; 
               }
               last_suffix = eptr->child; 
               if (jumps+edge_length(*eptr) <= 1) {
                  last_suffix->suffix_link = eptr->parent; 
                  last_suffix = NULL; 
               }
          
      
            }
         } else {
        
            add_leaf(*nptr, i, j, id);
            if (last_suffix != NULL) {
               last_suffix->suffix_link = nptr;
               last_suffix = NULL;  
            }      
         }
	   
         if (!END_EXPANSION_PHASE) {
            ++j;
            nptr = nptr->suffix_link; 
            if (jumps>0) { --jumps; }
         } 
         if (j > i) { 
            END_EXPANSION_PHASE = true; 
         }
      
    
      }
   }
} 




template <typename Sequence, 
          typename EndMarker, 

          typename EqualKey,
          typename Node,
          typename Edge> 

void 

_stree<Sequence, EndMarker, EqualKey, Node, Edge>
::clear(Node* nptr) {
   for (edge_iterator i = nptr->edges.begin(); i!=nptr->edges.end(); ++i) {
      clear(i->second->child);
      delete i->second;
   }
   delete nptr;
}




template <typename Sequence, 
          typename EndMarker, 

          typename EqualKey,
          typename Node,
          typename Edge> 

void

_stree<Sequence, EndMarker, EqualKey, Node, Edge>
::fix_label(Edge* eptr, SequenceID id) {
  
   if (eptr->sequence_id == id) { 
      edge_iterator i = eptr->child->edges.begin();
      Edge* low_eptr = i->second;
      eptr->sequence_id = low_eptr->sequence_id;
      eptr->first = low_eptr->first-edge_length(*eptr);
      eptr->last  = low_eptr->first-1;
   }
}




template <typename Sequence, 
          typename EndMarker, 

          typename EqualKey,
          typename Node,
          typename Edge> 

bool 

_stree<Sequence, EndMarker,  EqualKey, Node, Edge>
::erase_leaf(Edge* eptr, SequenceID id) {
   if (eptr->child->labels.find(id) != eptr->child->labels.end()) {
      if (eptr->child->labels.size() == 1) { 
         delete eptr->child;   
         delete eptr;         
         return true;
      } else {               
         eptr->child->labels.erase(id); 
         if (eptr->sequence_id == id) { 
            int elen = edge_length(*eptr);
            eptr->sequence_id = eptr->child->labels.begin()->first;
            eptr->last = lengths[eptr->sequence_id];
            eptr->first = eptr->last-elen+1;
         }
      }
   }
   return false;
}




template <typename Sequence, 
          typename EndMarker, 

          typename EqualKey,
          typename Node,
          typename Edge> 

void

_stree<Sequence, EndMarker,  EqualKey, Node, Edge>
::join_edge(Edge* eptr, SequenceID id) {
   Node* tmp = eptr->child;
  
   if (eptr->sequence_id == id) { 
      fix_label(eptr, id);
   }
  
   edge_iterator i = tmp->edges.begin();
   eptr->last  = eptr->last+edge_length(*(i->second));
   eptr->child = i->second->child;
  
   delete i->second;       
   delete tmp;      
}




template <typename Sequence, 
          typename EndMarker, 
          typename EqualKey,
          typename Node,
          typename Edge> 

void

_stree<Sequence, EndMarker,EqualKey, Node, Edge>
::erase(Node* nptr, SequenceID id) {
   std::vector<value_type> kill_list;
  
   for (edge_iterator i = nptr->edges.begin(); 
        i!=nptr->edges.end(); ++i) {
      if (!i->second->child->leaf) 
         erase(i->second->child, id);
    
      if (i->second->child->leaf) {
         if (erase_leaf(i->second, id)) {
            kill_list.push_back(i->first);
         }
      } else if (i->second->child->edges.size() == 1) {
         join_edge(i->second, id);
      } else {
         fix_label(i->second, id);
      }
   }
  
   for (typename std::vector<value_type>::iterator m = 
           kill_list.begin(); m != kill_list.end(); ++m) 
      nptr->edges.erase(*m);
 
   if (nptr->edges.size() == 0) {
      nptr->leaf = true;
      nptr->labels[id] = -1; 
   }
}

#endif
