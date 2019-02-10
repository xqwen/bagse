#include "classdef.h"
#include <fstream>
#include <sstream>



void show_banner(){

  fprintf(stderr, "\n\nBAGSE: BAyseian Gene Set Enriment Analysis\n\n");
  fprintf(stderr, "Usage:  bagse -d input_data_file -annot gene_set_annotation_file   [--load_zval]   [-lfdr lfdr_output_file]\n\n");
  
}



int main(int argc, char **argv){
  
  // creating the grid
  
  //olist.push_back(0.1);
  //phlist.push_back(0.05);

  char data_file[256];
  char gmap_file[256];
  char smap_file[256];
  char annot_file[256];
  char prior_dir[256];
  char output_pip[256];
  char output_wv[256];

  int csize = -1;
  int gsize = -1;
  int nthread = 1;
  int print_avg = 0;
  

  memset(data_file,0,256);


  char init_file[256];


  int find_egene = 0;
  int est = 0;
  double alpha = 0.05;

  double init_pi1 = 1e-3;

  double set_pi0 = -1;


  char ci_file[256];
  memset(ci_file,0,256); 
  memset(data_file,0,256); 
  memset(annot_file,0,256);
  memset(output_pip,0,256);

  int force_logistic = 0;
  int prob_annot = 0;
  
  int data_format = 1;
  
  double EM_thresh = 0.05;
  double l1_lambda = 0;
  double l2_lambda = 0;
 
  int use_ash = 1;

  for(int i=1;i<argc;i++){
    
    if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-data")==0){
      strcpy(data_file,argv[++i]);
      continue;
    }


    if(strcmp(argv[i], "-annot")==0 ){
      strcpy(annot_file,argv[++i]);
      continue;
    }
    
    if(strcmp(argv[i], "-t")==0 || strcmp(argv[i], "-thresh")==0){
      EM_thresh = atof(argv[++i]);
      continue;
    }
    
    if(strcmp(argv[i], "-l1_lambda")==0){
      l1_lambda = atof(argv[++i]);
      continue;
    }

    if(strcmp(argv[i], "-l2_lambda")==0){
      l2_lambda = atof(argv[++i]);
      continue;
    }
    

    
    if(strcmp(argv[i], "--load_bf")==0 || strcmp(argv[i], "--bf")==0 ){
      data_format = 2;
      continue;
    }
    
   
    if(strcmp(argv[i], "--load_zval")==0 || strcmp(argv[i], "--zval")==0){
      data_format = 3;
      continue;
    }
    
    if(strcmp(argv[i], "-lfdr")==0){
      strcpy(output_pip, argv[++i]);
      continue;
    }



    if(strcmp(argv[i], "-h")==0 || strcmp(argv[i], "-help")==0 ){
      show_banner();
      continue;
    }



    fprintf(stderr, "Error: undefined option %s\n", argv[i]);
    show_banner();
    exit(0);

  }    


  // checking mandatory arguments
  if(strlen(data_file)==0){
    fprintf(stderr,"Error: data file unspecified\n");
    show_banner();
    exit(0);
  }

    
  // a global variable 
  controller con;
  con.EM_thresh = EM_thresh;
  

  if(l1_lambda!=0){
    con.l1_lambda = l1_lambda;
    con.force_logistic = 1;
  }

  if(l2_lambda!=0){
    con.l2_lambda = l2_lambda;
    con.force_logistic =1;
  }
  
  
  con.use_ash = use_ash;

  switch(data_format){
  case 1:
    con.load_data(data_file);
    break;
  case 3:
    con.load_data_zscore(data_file);
    break;
  default:
    con.load_data(data_file);
    break;
  }
   

  con.load_annotation(annot_file);
  fprintf(stderr,"Initializing ... \n");
  con.estimate();
  if(strlen(output_pip)>0){
    //fprintf(stderr,"#-dump_pip output_pip file: %s \n",output_pip);
    con.dump_pip(output_pip);
  }

}


