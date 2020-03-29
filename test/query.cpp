#include "stree_query.h"


char* query_list;
char* query_dir;
char* ds_list;
char* ds_dir;

int bin_no = 2;
int lcs_min = 10;
int win_size = 3;
float epsilon = 0;


int openGapPenalty = -14;
int extensionGapPenalty  = -10;

int half_dm_size;
int dm_size;
int onehalf_dm_size;

using namespace std;

int main(int argc, char* argv[])
{
    extern char * optarg;
  int c;
   
  if(argc < 5 ) {
    cout << "usage: query \n";
    cout << "\t -l  <query proteins list> -REQUIRED \n";
    cout << "\t -d  <query proteins directroy> -REQUIRED \n";
    cout << "\t -L  <database proteins list> -REQUIRED \n";
    cout << "\t -D  <database proteins directory> -REQUIRED \n";
    cout << "\t -m  <the length threshold of maximal matches> \n";
    cout << "\t -e  <the distance threshold between symbols>  \n";
    cout << "\t -b  <the number of bins> \n";
    cout << "\t -o  <the opening gap penalty for dynamic programming> \n";
    cout << "\t -x  <the extension gap penalty for dynamic programming> \n";
    cout << "\t -w  <the size of window> \n";
    exit(1);
  }
  
  else{
    while ((c=getopt(argc,argv,"l:d:L:D:m:e:b:w:o:x:"))!=-1){
      switch(c){
      case 'l':
	query_list = optarg; 
	break;
      case 'd': 
	query_dir = optarg;
	break;
      case 'L': 
	ds_list = optarg;
	break;
      case 'D': 
	ds_dir = optarg;
	break;
      case 'm':
	lcs_min = atoi(optarg);
	break;
      case 'e':
	epsilon  = atof(optarg);
	break;
      case 'b':
	bin_no = atoi(optarg);
	break;
      case 'w':
	win_size = atoi(optarg);
	break;
      case 'o':
	openGapPenalty = atoi(optarg);
	break;
      case 'x':
	extensionGapPenalty = atoi(optarg);
	break;

      }
    }
  }
  if (query_list == ""){
    cout << "query list file name required\n";
    exit(1);
  }

  if (query_dir == ""){
    cout << "query dir required\n";
    exit(1);
  }

  if (ds_list == ""){
    cout << "database list file name required\n";
    exit(1);
  }

  if (ds_dir == ""){
    cout << "database dir required\n";
    exit(1);
  }

  half_dm_size = win_size - 1;
  dm_size = 2*half_dm_size;
  onehalf_dm_size = 3*half_dm_size;
  

  //check the argc
  //if(argc < 8){
  //cerr<<"Usage: ./query query_list_name query_pdb pdb_list_name pdb_dir bin_num distance lcs_min\n";
  //exit(EXIT_FAILURE);
  //}
  
  // output command line
  //cout <<endl << endl;
  //for(int i=0; i< argc; i++)
  //cout<< argv[i] << " ";
  // cout <<endl;

      
  //***********************************
  //win_size = 3;
  //half_dm_size = win_size -1;
  //***********************************
  
  
  //bin_no = atoi(argv[5]);
  //lcs_min = atoi(argv[7]);

  //float epsilon = atof(argv[6]);
  // if( bin_no*epsilon >= 1.0 )
  //epsilon = 1.0/((float)bin_no) - 0.0001;
  
  //MyRTreeQuery<seq_vector_of_int, vector_matches_list > q;
 
  stree_query<seq_vector_of_int, stree_matches_list > q;
  //q.query_stree(argv[1], argv[2],argv[3], argv[4], epsilon);
  q.query_stree(query_list, query_dir, ds_list, ds_dir, epsilon);
  

  //stree_query<seq_vector_of_int, seq_matches_list > q;
  //q.query(argv[1], argv[2],argv[3], argv[4], epsilon);


  //stree_query<seq_vector_vector_of_int, stree_matches_list > q;
  //q.query_stree(argv[1], argv[2],argv[3], argv[4], epsilon);
  
  //stree_query<seq_vector_vector_of_int, seq_matches_list > q; 
  //q.query(argv[1], argv[2],argv[3], argv[4], epsilon);


  //stree_query<seq_vector_of_string, stree_matches_list > q;
  //q.query_stree(argv[1], argv[2],argv[3], argv[4], epsilon);
  
  

  return 1;
}
