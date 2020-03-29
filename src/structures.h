
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <sstream>

#include "stree_map.h"

extern int bin_no;

extern int win_size;
extern int half_dm_size;


//#define win_size 3
//#define half_dm_size (win_size-1)

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

#ifndef __structure_h
#define __structure_h

using namespace std;

class Atom{
 public:

  //float	 coord[4];
  vector<float>  coord;
  //float dm[win_size][win_size];
  //float dmarray[half_dm_size];
  //vector<vector< float> > dm;
  vector<float>           dmarray;
  
  //int win_size;
  //int half_dm_size;
  
  
  Atom() {
    //for(int i=0; i<win_size ; i++)
    //dm[i][i]=0;
    coord.resize(3);
    //win_size = w;
    //half_dm_size = win_size - 1;
    //for(int i=0; i<win_size ; i++)
    //dm[i][i]=0;
    //dm.resize(win_size);
    
    // for(int i=0; i<win_size ; i++)
    //dm[i].resize(win_size);
    dmarray.resize(half_dm_size);

  }
  
  
  Atom(float x, float y, float z) {
    coord[0] = x;
    coord[1] = y;
    coord[2] = z;
    //for(int i=0; i<win_size ; i++)
    // dm[i][i]=0;
  }
  
  
  void set_coord(float x, float y, float z) {
    
    coord[0] = x;
    coord[1] = y;
    coord[2] = z;
  }
  
  
  
  void fill_dmarray(int pos, float distance){
    
    if(pos >= half_dm_size )
	  cerr << "distance array is out of bound index" << endl;
    dmarray[pos] = distance;
  }
  
    
  float pairwise_distance(Atom& neighbor){
    return sqrt( pow((coord[0]-neighbor.coord[0]),2) + pow((coord[1]-neighbor.coord[1]),2) + pow((coord[2]-neighbor.coord[2]),2) );
  }

  void print(){
    for (int i = 0; i < half_dm_size; i++)
      cout << dmarray[i] << " ";
    cout <<endl;
  }
    
} ;



class Amino_Acid
{
 public:
  
  Amino_Acid(){ 
    dmarray.resize(2*half_dm_size);
    anglearray.resize(half_dm_size);
  };
  
  
  int pid;
  int aid;
  
  Atom    CA_Atom;
  Atom    N_Atom;
  Atom    C_Atom;
  
  vector<float> dmarray;  
  vector<float> anglearray;
  
  vector<float> gammadata;
  vector<float> gammaarray;
  
  
  int arraysize() { 
    return dmarray.size(); 
    //return 2*half_dm_size; 
    //return 3*half_dm_size; 
  }
  
  /*
  void fill_dmarray(){
    
    int ai = 0;
    for(int i = 0; i < half_dm_size ; i++ ){
      dmarray[ai++] = CA_Atom.dmarray[i]/((win_size-1)*4.028);
      dmarray[ai++] = (anglearray[i]+1)/2;//cos = [-1,1]
    }
  }  
  */

  void fill_dmarray(){
    
    int ai = 0;
    for(int i = 0; i < half_dm_size ; ++i){

      float dist = CA_Atom.dmarray[i]/((win_size-1)*4.023);

      //cout <<  CA_Atom.dmarray[i] << " -> " << dist << endl;

      if(dist > 1.0)
	dmarray[ai++] = 0.999;
      else
	dmarray[ai++] = dist;
      


      float angle =  (anglearray[i]+1)/2;
      //cout <<  "-angle-" << anglearray[i] << " -> " << angle << endl;

      if(  angle > -100  && angle < 100 )
	dmarray[ai++] = angle;//cos = [-1,1]
      else
	dmarray[ai++] = 0.5;
    }
    
    int gsize =gammadata.size();
 
    //cout << "gsize---" << gsize <<endl;
    if( gsize > 0){
      for(int i=0; i < gsize; ++i){
	//cout << gammaarray[i] << " ";
        dmarray.push_back(gammadata[i]);
      }
      //cout <<endl; 
    }
  }  
  
  void fill_anglearray(int pos, float distance){
    if(pos >= half_dm_size)
      cerr << "angle array is out of bound index" << endl;
    anglearray[pos] = distance;
  }
  
  void set_aid(int id){ aid = id; }
  
  void direction(vector<float>& d)
    {
      d.resize(3);
      float x1=N_Atom.coord[0] - CA_Atom.coord[0];
      float y1=N_Atom.coord[1] - CA_Atom.coord[1];
      float z1=N_Atom.coord[2] - CA_Atom.coord[2];
     
      float x2=C_Atom.coord[0] - CA_Atom.coord[0];
      float y2=C_Atom.coord[1] - CA_Atom.coord[1];
      float z2=C_Atom.coord[2] - CA_Atom.coord[2];
      
      d[0] = y1*z2 - z1*y2;
      d[1] = -x1*z2 +z1*x2;
      d[2] = x1*y2 - x2*y1;
      
      float dist = sqrt(pow(d[0],2)+pow(d[1],2)+pow(d[2],2));
      d[0] = d[0]/dist;
      d[1] = d[1]/dist;
      d[2] = d[2]/dist;
      
    }
  
  
  
  void Compress(string& key) {	// object compression

      //char buff[10];
      //sprintf(&buff[0], "%d", pid);
      //memcpy(key, buff, sizeof(intb));
      //sprintf(&buff[0], "%d", aid);
      ///memcpy(key+sizeof(int), buff, sizeof(int));
      //memcpy(key+sizeof(int), &pid, sizeof(int));
      //memcpy(key+sizeof(int), &aid, sizeof(int));
      //memcpy(key, "1", 1);
      //memcpy(key+1, "2", 2);
      //cout << "--pid/aid-------" <<pid << "\t" << aid << endl;
      //memcpy(key+2*sizeof(int), "\0", 3);

      ostringstream os;
      os<< pid << " " << aid;
      key = os.str();
      // cout << "key is ----" << key << "+++" << endl;
  }

  friend ostream& operator<< (ostream& os, const Amino_Acid& aa)
    {
      os<< aa.pid << " " << aa.aid<< " ";
      //os<< atm.pdb_name << " " << atm.serial<< " ";
      //unsigned long i;
      //for (i = 0; i < half_dm_size; i++)
      //os << atm.dmarray[i] << " ";
      
      return os;
    }
    

} ;


class Structure{

 public:
  //string	    ident;
    int             pid;
    //int             win_size;
    //int             half_dm_size;
    
    vector<Amino_Acid>    vaa;
    
    Structure() {};
    
    Structure(int w) {
      //win_size = w;
      //half_dm_size = win_size - 1;
      //for(int i=0; i< vaa.size(); ++i)
      //vaa[i].set_win(w);
    };
    
    
    int size(){ return vaa.size() - win_size +1; };
    int arraysize(){ return 2*half_dm_size; };
    
    //Structure(char* id):  ident(id) {};
    //Structure(string id): ident(id) {};
    
    
    void assign_aid()
      {
	for(int i= 0; i<vaa.size(); i++){
	  vaa[i].pid = pid;
	  vaa[i].aid = i;
	}	
      }
    

    void print(){
      assign_aid();
      
      cout << "win_size " << win_size <<endl;
      
      for(int i= 0; i<  vaa.size() - win_size +1; ++i)
	{
	  
	  cout << "-----------atom coord------------ " <<endl;
	  cout<<vaa[i].aid << "\t" << vaa[i].CA_Atom.coord[0]  <<"\t" << vaa[i].CA_Atom.coord[1]  <<"\t" << vaa[i].CA_Atom.coord[2] <<endl;
  
	    //cout<< " -------distance matrix----- "<<endl;
	    //for(int row =0; row < win_size; row++)
	    //{
	    //for(int col = 0; col < win_size; col++)
	    //    cout<< vaa[i].CA_Atom.dm[row][col]<<"\t";
	    //cout<<endl;
	    //}

	  cout<< " -------distance array -----" <<endl;
	  //cout << "array size " << half_dm_size << endl;
	  
	  int asize = vaa[i].arraysize();
	  for(int rr =0; rr < asize; rr++)
	    cout<< vaa[i].dmarray[rr] <<"\t";
	  cout<<endl;
	  
	}
    }
    
} ;



//a specialized Sequence type
//combination of sequence and equal function

class seq_vector_of_int : public vector<int>
{

 private:
  int error;
  
 public:

  seq_vector_of_int() : error(0){};

  value_type endmarker() {return -1; }
  
  
  //-----------Euclidean distance-----------------
  ///*
  explicit seq_vector_of_int(float epsilon){
      //error = (int)(epsilon*bin_no);
      error = (int)(epsilon);
  }
  
  
  template<typename T>
  bool operator() (T m1, T m2)
    {
      if(error == 0 ) return m1==m2; // for exact matching
      // cout << ".";
      int d1, d2;
      div_t divresult1, divresult2;
      
      int sum =0 ;
      
      for(int i=0; i < 2*half_dm_size; i++){
	
	divresult1 = div (m1,bin_no);
	d1 = divresult1.rem;
	m1 = divresult1.quot;
	
	divresult2 = div (m2,bin_no);
	d2= divresult2.rem;;
	m2 = divresult2.quot;
	
	sum += (d1-d2)*(d1-d2);
	//if(abs(d1 - d2) > error)
	//return false;
	if(sum > error)
	  return false;
      }
      
      return sum <= error;
      
      //return true;
    }

  //*/



  //-----------Mahatan distance-----------------
  
  /*
  explicit seq_vector_of_int(float epsilon){
    error = (int)(epsilon * bin_no);
  }

  
  template<typename T>
  bool operator() (T m1, T m2)
    {
      if(error == 0 ) return m1==m2; // for exact matching
      // cout << ".";
      int d1, d2;
      div_t divresult1, divresult2;
      
      for(int i=0; i < 2*half_dm_size; i++){
	
	divresult1 = div (m1,bin_no);
	d1 = divresult1.rem;
	m1 = divresult1.quot;
	
	divresult2 = div (m2,bin_no);
	d2= divresult2.rem;;
	m2 = divresult2.quot;
	
	if(abs(d1 - d2) > error)
	  return false;
      }
      return true;
    }
  */
  
  
  
  template<typename Sequence>
  void Structure_to_Sequence(Structure& str, Sequence& seq){
    
    //int str_size = str.vaa.size() - win_size +1;
    int str_size = str.size();
    seq.resize(str_size); 
    
    for(int i = 0 ; i < str_size; ++i)
      seq[i] = mapping(str.vaa[i].dmarray); 
  }
  
  
  int mapping(vector<float>& dmarray){
    int mappedvalue = 0;
    int dsize = dmarray.size();
    for(int i=0; i < dsize; i++){
      mappedvalue += ((int)(dmarray[i]*bin_no))*((int)pow((float)bin_no, i));
      //cout << dmarray[i] << " ";
    }
    //cout << " -> " << mappedvalue << endl;
    
    return mappedvalue;
  }
  
  
  
  
  template<typename T>
  int similarity(T m1, T m2)
    {
      //if(error == 0 ) return m1==m2; // for exact matching
      // cout << ".";
      int d1, d2;
      div_t divresult1, divresult2;
      
      int sum =0 ;
      
      for(int i=0; i < 2*half_dm_size; i++){
	
	divresult1 = div (m1,bin_no);
	d1 = divresult1.rem;
	m1 = divresult1.quot;
	
	divresult2 = div (m2,bin_no);
	d2= divresult2.rem;;
	m2 = divresult2.quot;
	
	sum += (d1-d2)*(d1-d2);
	//if(abs(d1 - d2) > error)
	//return false;
	//if(sum > error)
	//return false;
      }

      return (int) (bin_no*((-1*sqrt((float)sum))+256*sqrt(2.0))/(256*sqrt(2.0)+40));
      
      //return true;
    }

  
};


//a specialized Sequence type
//not working because EndMarker have to be a constant.

class seq_vector_vector_of_int : public vector<vector<int> >
{
 protected:
  //typedef typename Sequence::value_type   value_type;  
 private:
  int error;
  
 public:
  seq_vector_vector_of_int() : error(0){};
  
  explicit seq_vector_vector_of_int(float epsilon){
      //error = (int)(epsilon * bin_no);
      error = (int)(epsilon);
  }

  value_type endmarker() {return vector<int>(2*half_dm_size, -1); }
  
  
  template<typename T>
  bool operator() (T m1, T m2)
    {
      if(error == 0 ) return m1 == m2;
      if(m1.size() != m2.size() ) {
	cout << "------different size-------" << m1.size() << " " << m2.size()<< endl;
	for(int i=0; i < m1.size() ; ++i)
	  cout << m1[i] <<" ";
	cout << endl;
	
	for(int i=0; i < m2.size() ; ++i)
	  cout << m2[i] <<" ";
	cout << endl;
	
	return false;
      }
      
      int size = m1.size();
      
      for(int i=0; i < size; ++i){
	if(abs(m1[i] - m2[i]) > error)
	  return false;
      }
      return true;
    }
  
  
  template<typename Sequence>
    void Structure_to_Sequence(Structure& str, Sequence& seq){
    
    int str_size = str.size() ;
    seq.resize(str_size);

    int arraysize = str.vaa[0].dmarray.size();
    //cout << arraysize << " ";
    
    for(int i = 0 ; i < str_size; ++i)
      //seq[i].resize(2*half_dm_size);
      seq[i].resize(arraysize);
    
    for(int i = 0 ; i < str_size; ++i)
      for(int j=0; j < arraysize; ++j)
	seq[i][j] = (int) (str.vaa[i].dmarray[j]*bin_no);
  }

  
  template<typename Sequence>
    void print(Sequence& seq) {
    
    int seq_size = seq.size();
    for(int i = 0 ; i < seq_size; ++i){
      int ssize = seq[i].size();  
      cout << "(";
      
    for(int j=0; j < ssize; ++j)
      cout << seq[i][j]<<  ",";
    cout << ")";
    }
    cout << endl;
    
  };  
  
};






class seq_vector_of_string : public vector<string>
{
 protected:
  //typedef typename Sequence::value_type   value_type;  
 private:
  int error;
  
 public:
  seq_vector_of_string() : error(0){};
  
  explicit seq_vector_of_string(float epsilon){
      //error = (int)(epsilon * bin_no);
      error = (int)(epsilon);
  }

  value_type endmarker() {return string(2*half_dm_size, '$'); }
  
  
  template<typename T>
  bool operator() (T m1, T m2)
    {
      if(error == 0 ) return m1 == m2;
      if(m1.size() != m2.size() ) {
	cout << "------different size-------" << m1.size() << " " << m2.size()<< endl;
	cout << m1 << "\t" << m2 <<endl;
	return false;
      }
      
      int size = m1.size();
      
      for(int i=0; i < size; ++i){
	if(abs(int(m1[i]) - int(m2[i])) > error)
	  return false;
      }
      return true;
    }
  
  
  template<typename Sequence>
    void Structure_to_Sequence(Structure& str, Sequence& seq){
    
    int str_size = str.size();
    seq.resize(str_size);
    
    int arraysize = str.vaa[0].dmarray.size();
    //cout << arraysize << " ";
    
    for(int i = 0 ; i < str_size; ++i)
      //seq[i].resize(2*half_dm_size);
      seq[i].resize(arraysize);
    
    for(int i = 0 ; i < str_size; ++i)
      for(int j=0; j < arraysize; ++j)
	seq[i][j] = char( 97+ int (str.vaa[i].dmarray[j]*bin_no));
  }

  
  template<typename Sequence>
    void print(Sequence& seq) {
    
    int seq_size = seq.size();

    for(int i = 0 ; i < seq_size; ++i){
    //int ssize = seq[i].size();  
      
      cout << "(";
      cout << seq[i];
      //for(int j=0; j < ssize; ++j)
    //cout << seq[i][j]<<  ",";
      cout << ")";
    }
    cout << endl;
    
  };  
  
};














// for seq_vector_vector_of_int
/*

*/

#endif

