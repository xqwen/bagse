#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "bagse.h"

void show_banner(){

    fprintf(stderr, "\n\nBAGSE: BAyseian Gene Set Enriment Analysis\n\n");
    fprintf(stderr, "Usage:  bagse -d input_data_file   [--load_zval]   [-fdr_out fdr_output_file]  [-fdr_level alpha]\n\n");

}



int main(int argc, char **argv){

    // creating the grid

    //olist.push_back(0.1);
    //phlist.push_back(0.05);

    char data_file[256];
    char output_fdr[256];
    int nthread = 1;
    int use_zval = 0;
    double EM_thresh = 0.1;


    memset(data_file,0,256);
    memset(output_fdr,0, 256);


    for(int i=1;i<argc;i++){

        if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-data")==0){
            strcpy(data_file,argv[++i]);
            continue;
        }


        if(strcmp(argv[i], "-t")==0 || strcmp(argv[i], "-thresh")==0){
            EM_thresh = atof(argv[++i]);
            continue;
        }


        if(strcmp(argv[i], "--load_zval")==0 || strcmp(argv[i], "--zval")==0){
            use_zval = 1;
            continue;
        }

        if(strcmp(argv[i], "-fdr_out")==0){
            strcpy(output_fdr, argv[++i]);
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
    BAGSE bagse;
    bagse.load_data(data_file, use_zval);
    bagse.run(EM_thresh);


}


