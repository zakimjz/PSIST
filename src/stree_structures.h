
#ifndef __STREE__STRUCTURES__
#define __STREE__STRUCTURES__

#include <vector>
#include <iostream>

extern int bin_no;

using namespace std;

class _matches 
{
 public:
  _matches(): lable(0), length(0), sequenceid(-1), score(0) {};
  
  float to_score()  {
    score = (float) length;
    return score;
  }
  
  long sequenceid;
  int lable;
  int length;
  float score;  
  int querylable;
  long queryid;
  
  void set_match(long sid,  int lab, int len, float s, int qlab)
    {
      sequenceid = sid;
      lable = lab;
      length = len;
      score = s ;
      querylable = qlab;
    }
};



//-------------------------------------------------------------
//-------------------------------------------------------------


class vector_matches_list
{
 public:
  vector<_matches> mlist;
  
  void clear() { mlist.clear(); }
  
  void add_matches(vector<_matches>& m){
    
    int m_size(m.size());
    //cout << "m.size--" << m_size <<endl;
    for(int i=0; i< m_size; i++){
      
      float score = m[i].to_score();
      long sid = m[i].sequenceid;
      if(sid + 1 > mlist.size() )
	mlist.resize(sid+1);
      
      if(score > mlist[sid].score){
	mlist[sid].score = score;
	mlist[sid].sequenceid = sid;
	//mlist[sid].set_match(sid,m[i].lable, m[i].length, score, querylable);
      }
      
    }
  }// end function


  
  void add_matches(vector<_matches>& m, int querylable){
    
    int m_size(m.size());
    //cout << "m.size--" << m_size <<endl;
    for(int i=0; i< m_size; i++){
      
      float score = m[i].to_score();
      long sid = m[i].sequenceid;
      if(sid + 1 > mlist.size() )
	mlist.resize(sid+1);
      
      if(score > mlist[sid].score){
	mlist[sid].score = score;
	mlist[sid].sequenceid = sid;
	//mlist[sid].set_match(sid,m[i].lable, m[i].length, score, querylable);
      }
      
    }
  }// end function
  
  void to_rank( vector<pair<int, float> >& r){
    
    int m_size = mlist.size();
    r.resize(m_size);
    for(int i=0; i< m_size; i++){
      r[i].first = mlist[i].sequenceid;
      r[i].second = mlist[i].score;
      //cout << "--" << r[i].first << "--" << r[i].second << endl;
    }
    clear();
  }

};



//-------------------------------------------------------------
//-------------------------------------------------------------

class lablelist
{
 public:

  lablelist(int l, int len, int ql) : lable(l),length(len),querylable(ql){};
  
  lablelist(const lablelist& l){
    lable = l.lable;
    length = l.length;
    querylable = l.querylable;
  }
  
  int lable;
  int length;
  int querylable;
  //int id;
};


class _lablelist_less
{
 public:
  bool operator() (const lablelist& a, const lablelist& b)
    {
      return a.length > b.length ;
    }
};



class slable
{
 public:
  int endpoint;
  int lableid;
  bool left;
  int hk;
  int lk;
  int length;
  
};

class _slable_less
{
 public:
  bool operator() (const slable& a, const slable& b)
    {
      return a.endpoint < b.endpoint ;
    }
};


class rectangle
{
 public:
  rectangle(int x1, int y1, int x2, int y2) : topx(x1), topy(y1), lowx(x2), lowy(y2){};
  
  int topx;
  int topy;
  int lowx;
  int lowy;
};



//-------------------------------------------------------------
//-------------------------------------------------------------

class seq_matches
{
 public:
  seq_matches() : sequenceid(-1),score(0){};
  long sequenceid; 
  float score;  
  vector<lablelist> ll;
  
  /*
  float to_score()
  {
    int ll_size = ll.size();
    int max_length = 0;
    for(int i=0; i< ll_size; ++i){
    if(max_length < ll[i].length)
	max_length = ll[i].length;
    }
    
    score = (float) max_length;
    return score;
  }
  */

  bool overlap(rectangle& query, rectangle& ref){
    if(query.lowx <= ref.topx && query.lowy <= ref.topy ) return false;
    if(query.topx >= ref.lowx && query.topy >= ref.lowy ) return false;
    return true;
  }
  

  float to_score(){
    
    int ll_size = ll.size();
    int max_length = 0;
    
    // vector<lablelist> ll2;
    //copy(ll.begin(), ll.end(), ll2.begin() );
    
    sort(ll.begin(), ll.end(), _lablelist_less() );
    vector<lablelist>::iterator iter;
    vector<lablelist>::iterator next;
    
    int gappenalty =  bin_no;
    //int gappenalty =  0;
    int gap =0;
    
    while(ll.size() > 0) {
      
      iter = ll.begin();
      max_length += iter->length;
      ++gap;
      rectangle r1(iter->querylable, iter->lable, iter->querylable + iter->length, iter->lable + iter->length );
      
      while(iter != ll.end() ){
	rectangle r2(iter->querylable, iter->lable, iter->querylable + iter->length, iter->lable + iter->length );
	
	if( overlap(r1, r2) == true ) {
	  next = ll.erase(iter);
	  iter = next;
	}
	else
	  ++iter;
      }
    }
    
    score = (float) (max_length - (gap-1)*gappenalty);
    return score;
    
  }
  

  






  /*
  float to_score(){
    
    //Horizational: query;  Vertical reference.
    int ll_size = ll.size();
    int max_length = 0;
    vector<slable> IL(2*ll_size); // <end point pos, rectangle id>
    vector<int>    V(ll_size, 0);
   
    map<int, pair<int, int> > mll;  
    
    for(int i=0; i< ll_size; ++i){
      IL[2*i].endpoint = ll[i].querylable;
      IL[2*i].lableid = i;
      IL[2*i].left    = true;
      IL[2*i].hk      = ll[i].lable;
      
      IL[2*i+1].endpoint = ll[i].querylable + ll[i].length - 1;
      IL[2*i+1].lableid = i;
      IL[2*i+1].left    = false;
      IL[2*i+1].lk      = ll[i].lable + ll[i].length - 1 ;
      
      mll[ IL[2*i+1].lk ] = make_pair(0,i);
      
    }

    sort(IL.begin(), IL.end(), _slable_less() );
    // cout << "IL size " << IL.size() << endl;
    
   
    map<int, pair<int, int> >::iterator iter;
    
    for(int i=0; i< 2*ll_size; ++i){
      
      if(IL[i].left == true){
	int hk = IL[i].hk;
	if( (iter = mll.upper_bound(hk) ) != mll.end() ){
	  V[IL[i].lableid] =  ll[IL[i].lableid].length + V[iter->second.second];
	  
	  //mll[IL[i].lableid].first = ll[IL[i].lableid].length + iter->second.first;
	}
      }
      else{
	
	int lk = IL[i].lk;
	int vk = V[IL[i].lableid];
	if(  (iter = mll.lower_bound(lk)) != mll.end() ){

	  if( (vk > V[iter->second.second]) && (iter->first <= lk) )
	    mll[lk] = make_pair(vk, IL[i].lableid);
	}
	for( iter = mll.begin(); iter != mll.end(); ++iter ){
	  if( (iter->first <= lk) && (V[iter->second.second] < vk ) )
	    mll.erase(iter);
	}
      }
    }
    
    
    //map<int, pair<int, int> >::iterator
    for(iter = mll.begin(); iter != mll.end(); ++iter){
      
      
      //iter = mll.end(); 
      //--iter;
      //iter = mll.begin();
      if( max_length < iter->second.first)
	max_length = iter->second.first;
      
    }
    
    score = (float) max_length;
    return score;
  }
  */  
  
  
};



class seq_matches_list 
{
 public:
  vector<seq_matches> mlist;
    
  void clear() { 
    mlist.clear();
    mlist.resize(0);
  }
  
  void add_matches(vector<_matches>& m){
    
    int m_size(m.size());
    for(int i=0; i< m_size; i++){
      
      float score = m[i].to_score();
      long sid = m[i].sequenceid;
      if(sid + 1 > mlist.size() )
	mlist.resize(sid+1);

      mlist[sid].sequenceid = sid;
      mlist[sid].ll.push_back(lablelist(m[i].lable, m[i].length, m[i].querylable));
    }
    
  }// end function
  
  
  void add_matches(vector<_matches>& m, int querylable){
    
    int m_size(m.size());
    for(int i=0; i< m_size; i++){
      
      float score = m[i].to_score();
      long sid = m[i].sequenceid;
      if(sid + 1 > mlist.size() )
	mlist.resize(sid+1);

      mlist[sid].sequenceid = sid;
      mlist[sid].ll.push_back(lablelist(m[i].lable, m[i].length, querylable));
    }
    
  }// end function
  
  
  void to_rank(vector<pair<int, float> >& r, int pos){
    
    int m_size = mlist.size();
    r.resize(m_size);
    for(int i=0; i< m_size; i++){
      r[i].first = mlist[i].sequenceid;
      r[i].second = mlist[i].to_score();
      //cout << "--" << r[i].first << "--" << r[i].second << endl;
    }
    clear();
  }
  
  
  void to_rank(vector<pair<int, float> >& r){
    
    int m_size = mlist.size();
    r.resize(m_size);
    for(int i=0; i< m_size; i++){
      r[i].first = mlist[i].sequenceid;
      r[i].second = mlist[i].to_score();
      //cout << "--" << r[i].first << "--" << r[i].second << endl;
    }
    clear();
  }

  void print(){
    int m_size = mlist.size();
    for(int i=0; i< m_size; ++i){
      cout << "--sid--" << mlist[i].sequenceid << "--score--" << mlist[i].to_score() << endl;
      int l_size = mlist[i].ll.size();
      for(int j=0; j<l_size ; ++j)
	cout << "(" << mlist[i].ll[j].lable << "," << mlist[i].ll[j].length << "," << mlist[i].ll[j].querylable << ")"  ;
      cout << endl;
    }
  }
  
  
};


//-------------------------------------------------------------
//-------------------------------------------------------------

class stree_matches_list 
{
 public:
  vector<vector<seq_matches> >  qmlist;
  
  void clear() { 
    qmlist.clear();
    qmlist.resize(0);
  }
    
  void add_matches(vector<_matches>& m){
    
    int m_size(m.size());
    //cout << "m.size--" << m_size <<endl;
    for(int i=0; i< m_size; i++){
      
      float score = m[i].to_score();
      long qid = m[i].queryid;
      long sid = m[i].sequenceid;
      if(qid + 1 > qmlist.size())
	qmlist.resize(qid+1);
      
      if(sid + 1 > qmlist[qid].size() )
	qmlist[qid].resize(sid+1);
      
      qmlist[qid][sid].sequenceid = sid;
      qmlist[qid][sid].ll.push_back(lablelist(m[i].lable, m[i].length, m[i].querylable));
    }
    
  }// end function

  
  void to_rank(vector<pair<int, float> >& r, int pos=0){
    
    if(pos >= qmlist.size() ) return;
    int m_size = qmlist[pos].size();
    r.resize(m_size);
    for(int i=0; i< m_size; i++){
      r[i].first = qmlist[pos][i].sequenceid;
      r[i].second = qmlist[pos][i].to_score();
      //cout << "--" << r[i].first << "--" << r[i].second << endl;
    }
    //clear();
  }

  
  void print(){
    int qm_size = qmlist.size();
    for(int i=0; i< qm_size; ++i){
      int m_size = qmlist[i].size();
      cout << "----query-----" << i << endl;
      
      for(int j=0; j< m_size; ++j){
      cout << "--sid--" << qmlist[i][j].sequenceid << "--score--" << qmlist[i][j].score << endl;
      int l_size = qmlist[i][j].ll.size();
      for(int k=0; k<l_size ; ++k)
	cout << "(" << qmlist[i][j].ll[k].lable << "," << qmlist[i][j].ll[k].length << "," << qmlist[i][i].ll[k].querylable << ")"  ;
      cout << endl;
      }
    }

  }
  
  
};


//-------------------------------------------------------------
//-------------------------------------------------------------


//template <typename T>
/*
class _matches_less
{
 protected:
  //typedef T match_type;
  typedef _matches match_type;
 public:
  bool operator() (const match_type& a, const match_type& b)
    {
      return a.score > b.score ;
    }
};
*/


#endif
