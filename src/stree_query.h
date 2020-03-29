
#include <stdlib.h>
#include "manipulate_structures.h"
#include "qrank.h"
#include "TimeTracker.h"

#define batch_num 1

#define _DP_


//extern int bin_no;
extern int lcs_min;

//extern int win_size;
//extern int half_dm_size;

//extern int win_size;

using namespace std;

#ifndef __STREE_QUERY_H
#define __STREE_QUERY_H


template <typename Sequence, typename match_list_type> 
class stree_query
{
 protected:
  typedef typename Sequence::difference_type difference_type;
  typedef typename Sequence::size_type   size_type;
  typedef typename Sequence::value_type   value_type;
  typedef long SequenceID;
  
  typedef _matches match_type;
  typedef Sequence seq_type;
  
  //typedef stree_map<Sequence, -1> tree_type;
  //typedef stree_map<Sequence, -1, seq_type> inexact_tree_type;
  
  typedef stree_map<Sequence, value_type> tree_type;
  typedef stree_map<Sequence, value_type, seq_type> inexact_tree_type;
  
  
  //typedef  Matches match_list_type;
  
  
 public:
  
  qrank  qrk, dprank;
  match_list_type _m;
  
  string stree_name;
  
  stree_query() {};
  stree_query(char* stree_n) : stree_name(stree_n) {};
  ~stree_query(){};
  
  void insert(tree_type* t, string pdb_name, long sid){
    
    //string gammaname = pdb_name + ".gamma";
    pdb_name += ".coord";
    
    Structure s = read_parsed_pdb_to_structure(pdb_name.c_str());
    
    if(s.vaa.size() > 0 ){
      
      //---add seq info---------
      //read_gamma_to_structure(gammaname.c_str(), s);
      
      calculate_distance_matrix(s);
      calculate_dm_array(s);

      //---add seq info---------
      //calculate_gamma_array(s);
      
      //cout << pdb_name << endl;
      
      //s.print();
      Sequence seq;
      seq.Structure_to_Sequence(s, seq);
      //cout << "****************************s" << endl;
      
      t->insert(seq, sid);
    }
    else
      cout << "------empty pdb-----------" << pdb_name << endl;
    
  }
  
  
  void build_stree(tree_type* tree, char* pdblist, char *pdb_dir)
    {  
      ifstream fin(pdblist, ios::in);
      if(!fin){
	cerr << "couldnot open list "<< endl;
	exit(0);
      }
      if( opendir(pdb_dir) == NULL) {
	fprintf(stderr,"cannot open directory: %s\n", pdb_dir);
	return;
      }
      
      //tree = new stree_map<Sequence, -1>();
      long sid =0;
      
      while (fin)
	{  
	  if (! fin.good()) continue; // skip newlines, etc.
	  
	  pdbscop p;
	  fin >> p.name >> p.cl >> p.fold >> p.sf >> p.fam;
	  if(p.name.size()  < 1 ) continue;
	  
	  //qrk.add_pdb(p);
	  
	  //string pdb_name = pdb_dir + p.name + ".coord";
	  string pdb_name = pdb_dir + p.name;
	  insert(tree, pdb_name, sid);
	  sid++;
	} 	
      fin.close() ;
    }
  
  
  void dynprog(Sequence q, inexact_tree_type* tree, vector<pair<int, float> >& rank, float epsilon){
    
    int out_num = topk;
    //int out_num = knn;
    
    rank.resize(out_num);
    //int rsize = rank.size();
    //for(int i=0; i< rsize; ++i){
    for(int i=0; i< out_num; ++i){
      if(rank[i].first >=0 )
	rank[i].second = (float) localAlign(q, tree->sequences[rank[i].first], epsilon);
    }
  }
  

  void  dynprog(Structure& q, qrank& qrk, char *pdb_dir, vector<pair<int, float> >& rank){
  
    int out_num = topk;
    //int out_num = knn;
    rank.resize(out_num);
    
    for(int i=0; i< out_num; ++i)
      if(rank[i].first >=0 ){
	
	string pdb_name = pdb_dir + qrk.pdbs[rank[i].first].name;
	Structure s = read_parsed_pdb_to_structure(pdb_name.c_str());
	if(s.vaa.size() >0 ){
	  calculate_distance_matrix(s);
	  calculate_dm_array(s);
	  rank[i].second = localAlign(q,s);
	}
	else
	  continue;
      }
    
  }



  void  dynprog(string qname, qrank& qrk, char *pdb_dir, vector<pair<int, float> >& rank){
    
    int out_num = topk;
    //int out_num = knn;
    rank.resize(out_num);
    
    //cout <<  " read ---------" << qname << endl;
    
    Structure q = read_parsed_pdb_to_structure(qname.c_str());
    //cout << q.size() << endl;
    //cout << "--queyr- " << qname << "-size-"  << q.size() << endl;
    if(q.vaa.size() >0 ){
      calculate_distance_matrix(q);
      calculate_dm_array(q);
    }
    else return;
    
    for(int i=0; i< out_num; ++i)
      if(rank[i].first >=0 ){
	
	string pdb_name = pdb_dir + qrk.pdbs[rank[i].first].name+ ".coord";
	Structure s = read_parsed_pdb_to_structure(pdb_name.c_str());
	//cout << "--db- " << pdb_name << "-size-"  << s.size() << endl;
	
	if(s.vaa.size() >0 ){
	  calculate_distance_matrix(s);
	  calculate_dm_array(s);
	  rank[i].second = localAlign(q,s);
	  //cout << "--end dyn--" <<endl;
	}
	else
	  continue;

      }
    
  }
  
  
  
  /*
  void query_stree(char* query_list, char* query_dir, char* pdb_list, char* pdb_dir, float epsilon){
    
    TimeTracker timeIndex, timeQuery;
    
    vector<int> wins;
    wins.push_back(3);
    //wins.push_back(5);
    //wins.push_back(7);

    vector<tree_type*> trees;
    
    timeIndex.start();
    for(int i=0; i< wins.size(); ++i){
      tree_type* ptree = new tree_type();
      win_size = wins[i];
      half_dm_size = win_size -1;
      build_stree(ptree, pdb_list, pdb_dir);
      trees.push_back(ptree);
    }
    
    timeIndex.stop();


    cout << "-------------end construction------------" <<endl;
    


    timeQuery.start();
    qrk.read_pdbs(pdb_list);
    dprank.read_pdbs(pdb_list);
    
    //build_stree(qtree, query_list, query_dir);
    
    vector<inexact_tree_type*> rtrees;
    for(int i=0; i< wins.size(); ++i){
      inexact_tree_type *rptree;
      rptree= (inexact_tree_type* ) trees[i];
      rtrees.push_back(rptree);
    }
    
    ifstream fin(query_list, ios::in);
    if(!fin){
      cerr << "couldnot open list "<< endl;
      exit(0);
    }
    if( opendir(query_dir) == NULL) {
      fprintf(stderr,"cannot open directory: %s\n", pdb_dir);
      return;
    }

    //tree = new stree_map<Sequence, -1>();
    long sid =0;
    int count = 0;
    vector<pdbscop> vp;
   
    while (fin){  
      if (! fin.good()) continue; // skip newlines, etc.
      
      pdbscop p;
      fin >> p.name >> p.cl >> p.fold >> p.sf >> p.fam;
      if(p.name.size()  < 1 ) continue;
      
      // vp.push_back(p);
      
      cout << sid+1 << "--query--" <<p.name << endl;      
      
      //string pdb_name = query_dir + p.name + ".coord";
      string pdb_name = query_dir + p.name;
 
      vector<tree_type*> qtrees;
      for(int i=0; i< wins.size(); ++i){
	tree_type* qtree = new tree_type();
	qtrees.push_back(qtree);
      }
      
      for(int i=0; i< wins.size(); ++i){
	win_size = wins[i];      
	half_dm_size = win_size -1;
	insert(qtrees[i], pdb_name, count);
      }
      
      ++sid;
      //++count;
      
      
 
      // if(count == batch_num){
      
      inexact_tree_type *rqtree;
      //rqtree= (inexact_tree_type* ) qtree;

      vector<pair<int, float> > ranksum;

      
      for(int w=0; w< wins.size(); ++w){
	
	win_size = wins[w];
	half_dm_size = win_size -1;
      
	vector<match_type> res;
	rqtree= (inexact_tree_type* ) qtrees[w];
	

	//if(epsilon > 0)
	rtrees[w]->find_appr_lcs_stree(rqtree->root,rqtree->sequences, lcs_min, res, epsilon);
	//else
	//ptree->find_lcs_stree(qtree->root,qtree->sequences, lcs_min, res);
	
	_m.add_matches(res);
	
	//for(int i=0; i< count; ++i){
	vector<pair<int, float> > r;
	_m.to_rank(r, 0);
	

	if(w == 0){
	  ranksum.resize(r.size());
	  for(int k=0; k< r.size(); ++k){
	    ranksum[k].first = k;
	    ranksum[k].second = r[k].second;
	  }
	}
	else{
	  for(int k=0; k< r.size(); ++k)
	    ranksum[k].second += r[k].second;
	}
	
      }
      
      qrk.add_to_perct_hits(p, ranksum);
      //qrk.print_arank(vp[i], r);
	
      win_size = wins[0];
      half_dm_size = win_size -1;
      //dynprog(qtrees[0]->sequences[0], rtrees[0], ranksum, epsilon);
      //dprank.add_to_perct_hits(p, ranksum);
      //qrk.print_arank(vp[i], r);
     
	
      //_m.print();
      //vp.clear();
      //vp.resize(0);
      _m.clear();
      count = 0;
      //delete qtrees;
      //qtree = new tree_type();
      qtrees.clear();
      qtrees.resize(wins.size());
    }
    
    
    
    
    fin.close();
    timeQuery.stop();
    
    //qrk.set_query_num(query_num);
    qrk.print_perct_hits();
    cout << "-------------dynprog---------------------------------" <<endl;
    //dprank.print_perct_hits();
    
    cout << endl;
    cout << "----------------------------------------------" <<endl;
    cout << "--index time--" << timeIndex.getTotalTime() << endl;
    cout << "--query time--" << timeQuery.getTotalTime() << "--Avg--" << timeQuery.getTotalTime()/((float)qrk.get_query_num() ) << endl;
    
  }
  

  
  */  

  
  ///*
  
  void query_stree(char* query_list, char* query_dir, char* pdb_list, char* pdb_dir, float epsilon){

      //win_size = 3;
      //half_dm_size = win_size -1;
    
    TimeTracker timeIndex, timeQuery;
    
    timeIndex.start();
    //tree_type* ptree = new tree_type(-1);
    tree_type* ptree = new tree_type(Sequence().endmarker() );
    build_stree(ptree, pdb_list, pdb_dir);
    timeIndex.stop();

    timeQuery.start();
    qrk.read_pdbs(pdb_list);
    dprank.read_pdbs(pdb_list);
    
    //build_stree(qtree, query_list, query_dir);

    inexact_tree_type *rptree;
    rptree= (inexact_tree_type* ) ptree;
    
    ifstream fin(query_list, ios::in);
    if(!fin){
      cerr << "couldnot open list "<< endl;
      exit(0);
    }
    if( opendir(query_dir) == NULL) {
      fprintf(stderr,"cannot open directory: %s\n", pdb_dir);
      return;
    }

    //tree = new stree_map<Sequence, -1>();
    long sid =0;
    int count = 0;
    //int tree_num = 0;

    vector<pdbscop> vp;
    
    tree_type* qtree;

    //qtree = new tree_type(-1);
    qtree = new tree_type(Sequence().endmarker() );
    
    while (fin){  
      if (! fin.good()) continue; // skip newlines, etc.
      
      pdbscop p;
      fin >> p.name >> p.cl >> p.fold >> p.sf >> p.fam;
      if(p.name.size()  < 1 ) continue;
      
      vp.push_back(p);
      
      //cout << sid+1 << "--query--" <<p.name << endl;      
      
      //string pdb_name = query_dir + p.name + ".coord";
      string pdb_name = query_dir + p.name;

      insert(qtree, pdb_name, count);
      
      ++sid;
      ++count;
      
 
      if(count == batch_num){
	
  	inexact_tree_type *rqtree;
	rqtree= (inexact_tree_type* ) qtree;
	
	vector<match_type> res;
	//if(epsilon > 0)
	  rptree->find_appr_lcs_stree(rqtree->root,rqtree->sequences, lcs_min, res, epsilon);
	  //else
	  //ptree->find_lcs_stree(qtree->root,qtree->sequences, lcs_min, res);
	
	_m.add_matches(res);
	
	for(int i=0; i< count; ++i){
	  vector<pair<int, float> > r;
	  _m.to_rank(r, i);
	  
	  //int rsize = r.size();
	  //for(int k=0; k < rsize; ++k)
	  //if(r[k].first >= 0)
	  //  r[k].second = r[k].second/((float)rptree->sequences[r[k].first].size());
	  //else
	  //  r[k].second =0;
	    
	  
	  //sort(r.begin(), r.end(), _rank_less() );
	  qrk.add_to_perct_hits(vp[i], r);
	  //qrk.print_arank(vp[i], r);

#ifdef _DP_
	  dynprog(rqtree->sequences[i], rptree, r, epsilon);

	  //string query_name = query_dir + vp[i].name +".coord"; 
	  //dynprog(query_name, qrk, pdb_dir, r);

	  dprank.add_to_perct_hits(vp[i], r);
#endif
	  //qrk.print_arank(vp[i], r);
	  
	} // end for i
	
	
	//_m.print();
	vp.clear();
	vp.resize(0);
	_m.clear();
	count = 0;
	delete qtree;
	//qtree = new tree_type(-1);
	qtree = new tree_type(Sequence().endmarker());
	
      }
      

    }
    
    // remaining pdbs--------------------
    if(count >0){
      inexact_tree_type *rqtree;
      rqtree= (inexact_tree_type* ) qtree;
      
      vector<match_type> res;
      //if(epsilon > 0 )
	rptree->find_appr_lcs_stree(rqtree->root,rqtree->sequences, lcs_min, res, epsilon);
	//else
	//ptree->find_lcs_stree(qtree->root,qtree->sequences, lcs_min, res);
      
      _m.add_matches(res);
      
      for(int i=0; i< count; ++i){
	vector<pair<int, float> > r;
	_m.to_rank(r, i);

	//int rsize = r.size();
	//for(int k=0; k < rsize; ++k)
	//if(r[k].first >= 0)
	//  r[k].second = r[k].second/((float)rptree->sequences[r[k].first].size());
	//else
	//  r[k].second =0;
	
	qrk.add_to_perct_hits(vp[i], r);
	//qrk.print_arank(vp[i], r);

#ifdef _DP_
	dynprog(rqtree->sequences[i], rptree, r, epsilon);
	
	//string query_name = query_dir + vp[i].name+".coord";
	//dynprog(query_name, qrk, pdb_dir, r);
	  
	dprank.add_to_perct_hits(vp[i], r);	
	//qrk.print_arank(vp[i], r); 
#endif
      }
      //_m.print();
      
    }
    

    //--------------------------------------------------
    
    
    fin.close();
    timeQuery.stop();
    
    //qrk.set_query_num(query_num);
    qrk.print_perct_hits();

#ifdef _DP_
    cout << "-------------dynprog---------------------------------" <<endl;
    dprank.print_perct_hits();
#endif
    
    cout << endl;
    cout << "----------------------------------------------" <<endl;
    cout << "--index time--" << timeIndex.getTotalTime() << endl;
    cout << "--query time--" << timeQuery.getTotalTime() << "--Avg--" << timeQuery.getTotalTime()/((float)qrk.get_query_num() ) << endl;
  
  }

  //*/ 

  /*
  void query(char* query_list, char* query_dir, char* pdb_list, char* pdb_dir, float epsilon){

    win_size = 4;
    half_dm_size = win_size -1;
    
    TimeTracker ruIndex, ruQuery;
    
    ifstream fin(query_list, ios::in);
    if(!fin){
      cerr << "couldnot open querylist "<< endl;
      exit(0);
    }
    if( opendir(query_dir) == NULL) {
      fprintf(stderr,"cannot open directory: %s\n", query_dir);
      return;
    }
    
    ruIndex.start();
    //tree_type* ptree = new tree_type(-1);
    tree_type* ptree = new tree_type(Sequence().endmarker());
    
    build_stree(ptree, pdb_list, pdb_dir);
    ruIndex.stop();
    
    ruQuery.start();
    qrk.read_pdbs(pdb_list);
    dprank.read_pdbs(pdb_list);

    inexact_tree_type* tree1;
    
    //if(epsilon >0 )
    tree1 = (inexact_tree_type* ) ptree;
    tree1->emarker = Sequence().endmarker();
    
    
    int query_num =0;
    while (fin)
      {  
	if (! fin.good()) continue; // skip newlines, etc.
	
	pdbscop q;
	fin >> q.name >> q.cl >> q.fold >> q.sf >> q.fam;
	if(q.name.size()  < 1 ) continue;
	
	string query_name = query_dir + q.name + ".coord";   
	//string query_name = query_dir + q.name; 

	Structure s = read_parsed_pdb_to_structure(query_name.c_str());
	//cout << q.name << "\t" << s.vaa.size() << endl;
	cout << query_num << "--query--" <<q.name << endl;
	
	
	//***********************************
        //win_size = 3;
	//half_dm_size = win_size -1;
	//***********************************

	if(s.vaa.size() >= win_size ){
	  
	  calculate_distance_matrix(s);
	  calculate_dm_array(s);
	  //s.print();
	  
	  Sequence seq, seq2;
	  seq.Structure_to_Sequence(s, seq);
	  //copy(seq.begin(), seq.end(), seq2.begin());
	  seq2.Structure_to_Sequence(s, seq2);
	  
	  int querylable =0;
	  // value_type pref = -1;
	  value_type pref = seq.endmarker();
	
	
	  
	  while(seq.size() >= lcs_min){
	    
	    vector<match_type> res;
	    if(epsilon > 0)
	      tree1->find_appr_lcs(seq, lcs_min, res, epsilon, pref);
	    else
	      ptree->find_lcs(seq, lcs_min, res, pref);
	    
	
	    _m.add_matches(res, querylable);
	    pref = *(seq.begin());
	    seq.erase(seq.begin());
	    ++querylable;
	    // cout << "seq size " << seq.size() <<endl;
	  }
	  
	 
	  //_m.print();
	  vector<pair<int, float> > r;
	  _m.to_rank(r, 0);
	  //qrk.print_arank(q,r);
	  qrk.add_to_perct_hits(q, r);
	  
	  dynprog(seq2, tree1, r, epsilon);
	  //dynprog(s, qrk, pdb_dir, r);
	  
	  dprank.add_to_perct_hits(q, r);
	  
	  ++query_num;
	  
	}
	
      }//end while
 
    fin.close();
    ruQuery.stop();
    
    //qrk.set_query_num(query_num);
    qrk.print_perct_hits();
    cout << endl;
    cout << "----------------------------------------------" <<endl;
    dprank.print_perct_hits();
    
    cout << endl;
    cout << "----------------------------------------------" <<endl;
    cout << "--index time--" << ruIndex.getTotalTime() << endl;
    cout << "--query time--" << ruQuery.getTotalTime() << "--Avg--" << ruQuery.getTotalTime()/((float)query_num) << endl;
  
  }
  */
  void load_stree(){}
  
  
};




#endif
