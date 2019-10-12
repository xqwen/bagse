using namespace std;

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <gsl/gsl_cdf.h>
#include "bagse_logistic.h"


void BAGSE_logistic::load_data(char *filename, int use_zval, char *annot_file){



    ifstream dfile(filename);
    string line;
    istringstream ins;


    string loc_id;
    double beta;
    double se_beta;
    string annot;

    int loc_count = 0;

    vector<double> beta_vec;
    vector<double> se_vec;


    // for grid setting
    double min = 1; 
    double max = 0;; 
    while(getline(dfile,line)){

        ins.clear();
        ins.str(line);

        ins>>loc_id>>beta;

        if(!use_zval)
            ins>>se_beta;
        else
            se_beta = 1;

        if(!use_zval)
            ins>>se_beta;
        else
            se_beta = 1;

        if(use_zval == -1) {// pvalue is used
            if(beta<1e-16)
                beta = 1e-16;
            beta = gsl_cdf_ugaussian_Qinv (beta/2);
        }

        loc_vec.push_back(loc_id);
        gene_hash[loc_id] = 100;
        beta_vec.push_back(beta);
        se_vec.push_back(se_beta);

        
        if(se_beta <= min)
            min = se_beta;

        if(beta*beta - se_beta*se_beta > max){
            max = pow(beta,2)-pow(se_beta,2);
        }

    }

    dfile.close();

    // compute BF vector
    double phi_min = min/10;
    double phi_max = 2*sqrt(max);
    if(phi_max<phi_min){
        phi_max = 8*phi_min;
    }

    vector<double> grid_vec = make_grid(phi_min, phi_max);
    grid_size = grid_vec.size();



    N = beta_vec.size();


    fprintf(stderr, "Initializing ... \n");
    //fprintf(stderr, "N=%d\t K=%d\t L=%d\t  M=%d\n", N, grid_size, annot_size, K);

    // compute BF factors
    for(int i=0;i<N;i++){
        vector<double> bf_vec;
        for(int j = 0; j<grid_size+1;j++){

            double log10_bf = 0;
            if(j>0){
                log10_bf = compute_log10_BF(beta_vec[i], se_vec[i],grid_vec[j-1]);
            }
            //printf("index = %d log10_bf = %f %d %d\n", index, log10_bf, annot_vec[i], j);
            bf_vec.push_back(log10_bf);
        }


        log10_BF_matrix.push_back(bf_vec);  
    }
    /*
    for(int i=0;i<N;i++){
        for(int j=0;j<grid_size+1;j++){
            printf("%7.3f  ", log10_BF_matrix[i][j]);
        }
        printf("\n");
    }
    printf("N = %d grid_size = %d\n",N, grid_size);

    */



    load_annotation(annot_file);
}



void BAGSE_logistic::load_annotation(char *annot_file){

    map<string, vector<double> > annot_map;

    map<int, int> col2cat;
    map<int, int> col2cpos;
    map<int, int> col2dpos;
    int col_count = -1;



    ifstream afile(annot_file);
    string line;
    istringstream ins;


    // parse header

    getline(afile,line);
    ins.clear();
    ins.str(line);
    string token;
    
    
    kd = 0; // count of discrete annotation
    kc = 0; // count of continuous annotation
    Xc = 0;
    Xd = 0;
    dlevel = 0;


    while(ins>>token){

        if(col_count==-1){
            if(token == "Gene" || token == "gene" || token == "GENE"){
                col_count = 0;
                continue;
            }else{
                break;
            }

        }

        string cat = token.substr(token.size()-2, 2);
        // continuous
        if(cat == "_c" || cat =="_C"){
            col2cat[col_count] = 1;
            col2cpos[col_count] = kc; // position in Xc matrix
            string name = token.substr(0,token.size()-2);
            cvar_name_vec.push_back(name);
            kc++;
        } // discrete/categorical
        else{
            col2cat[col_count] = 2;
            col2dpos[col_count] = kd; // position in Xd matrix
            dvar_name_vec.push_back(token);
            kd++;
        }

        col_count++;
    }


    if(col_count == -1){
        fprintf(stderr,"\nError: invalid header in the annotation file.\n");
        exit(1);
    }

    // read in data

    while(getline(afile,line)){

        ins.clear();
        ins.str(line);
        string gene;


        ins>>gene;

        if(gene_hash[gene] != 100)
            continue;

        double val;
        vector<double> avec;
        while(ins>>val){
            avec.push_back(val);
        }
        annot_map[gene] = avec;
    }

    afile.close();

    if(kc>0){
        Xc = gsl_matrix_calloc(N,kc);
    }


    if(kd>0){
        Xd = gsl_matrix_int_calloc(N,kd);
        dlevel = gsl_vector_int_calloc(kd);
    }




    for (int i=0;i<N;i++){

        string gene_id = loc_vec[i];

        vector<double> avec((kd+kc),0.0);
        if(annot_map.find(gene_id)!=annot_map.end())
            avec = annot_map[gene_id];



        for(int k=0;k<avec.size();k++){

            if(col2cat[k] == 1){
                gsl_matrix_set(Xc,i,col2cpos[k],avec[k]);
            }else{
                gsl_matrix_int_set(Xd, i, col2dpos[k],int(avec[k]));
            }

        }

    }

    for(int i=0;i<kd;i++){
        gsl_vector_int_set(dlevel,i,count_factor_level(i));
    }

    

}



int BAGSE_logistic::count_factor_level(int col){

    map<int, int> rcd;
    for(int i=0;i<N;i++){
        int val = gsl_matrix_int_get(Xd,i,col);
        rcd[val] = 1;
    }

    return rcd.size();
}




vector<double> BAGSE_logistic::make_grid(double phi_min, double phi_max){

    vector<double> gvec;
    gvec.push_back(phi_max); //phi_max
    double phi = phi_max;
    while(phi>phi_min){
        phi = phi/sqrt(2);
        gvec.push_back(phi); //phi
    }
    std::sort(gvec.begin(),gvec.end());

    return gvec;
}




void BAGSE_logistic::run(double thresh){  

    // initial setting
    vector<double> wts_vec(grid_size, 1.0/grid_size );
    fprintf(stderr, "Running EM ... \n\n");

    int ncoef = 0;
    for(int i=0;i<kd;i++){
        ncoef += gsl_vector_int_get(dlevel,i)-1;
    } 

    ncoef += 1+kc;
    
    gsl_vector *alpha_init_vec= gsl_vector_calloc(ncoef);

    // run EM
    final_loglik = gem.EM_run(log10_BF_matrix, alpha_init_vec, wts_vec, Xd, dlevel, Xc, thresh);


    gsl_vector *alpha_test;
    double loglik;
    //report results
    gsl_vector *alpha_est = gem.get_alpha_estimate();
    double est = gsl_vector_get(alpha_est,0);
    printf("%15s  %7.3f   ","Intercept", est);

    alpha_test = gsl_vector_calloc(alpha_est->size);
    gsl_vector_memcpy(alpha_test, alpha_est);
    gsl_vector_set(alpha_test,0,0);
    loglik = gem.compute_loglik(alpha_test);
    gsl_vector_free(alpha_test);
    double se = sqrt(pow(est,2)/(2*(final_loglik-loglik)));
    printf("  %7.3f %7.3f\n", est-1.96*se, est+1.96*se);


    int index = 1;
    for(int i=0;i<dvar_name_vec.size();i++){
        int level = gsl_vector_int_get(dlevel,i);
        string prefix = dvar_name_vec[i];
        for(int j=1;j<level;j++){
            ostringstream stream;
            stream <<prefix<<"."<<j;
            string label = stream.str();
            est = gsl_vector_get(alpha_est, index);
            printf("%15s  %7.3f   ",label.c_str(),est);

            alpha_test = gsl_vector_calloc(alpha_est->size);
            gsl_vector_memcpy(alpha_test, alpha_est);
            gsl_vector_set(alpha_test,index,0);
            loglik = gem.compute_loglik(alpha_test);
            gsl_vector_free(alpha_test);
            se = sqrt(pow(est,2)/(2*(final_loglik-loglik)));
            printf("  %7.3f %7.3f\n", est-1.96*se, est+1.96*se);
            index++;
        }
    }



    for(int i=0;i<cvar_name_vec.size();i++){
        string label =  cvar_name_vec[i];
        double est = gsl_vector_get(alpha_est, index);
        printf("%15s  %7.3f   ",label.c_str(),est);

        alpha_test = gsl_vector_calloc(alpha_est->size);
        gsl_vector_memcpy(alpha_test, alpha_est);
        gsl_vector_set(alpha_test,index,0);
        loglik = gem.compute_loglik(alpha_test);
        gsl_vector_free(alpha_test);
        se = sqrt(pow(est,2)/(2*(final_loglik-loglik)));
        printf("  %7.3f %7.3f\n", est-1.96*se, est+1.96*se);
        index++;
    }


    gsl_vector_free(alpha_init_vec);
    fprintf(stderr,"\n\n");

} 



void BAGSE_logistic::fdr_control(char *fdr_out, double fdr_level){

    gsl_vector * pip_vec = gem.get_pip_vec();
    vector<double> lfdr_vec;
    for(int i=0;i<N;i++){
        lfdr_vec.push_back(1 - gsl_vector_get(pip_vec,i));
    }

    vector<double> lfdr_sort_vec = lfdr_vec;
    std::sort(lfdr_sort_vec.begin(), lfdr_sort_vec.end());
    double cutoff = 0;
    double sum = 0;
    for(int i=0;i<N;i++){
        sum += lfdr_sort_vec[i];
        if(sum/(i+1) > fdr_level){
            fprintf(stderr, "\n\n   FDR control: %d genes rejected at level %.2f\n\n", i, fdr_level);
            break;
        }
        cutoff = lfdr_sort_vec[i];
    }

    FILE *fout;
    if(strlen(fdr_out) == 0){
        fout = fopen("fdr.out", "w");
    }else
        fout = fopen(fdr_out,"w");

    for(int i=0;i<N;i++){
        int rej = 0;
        if(lfdr_vec[i]<=cutoff){
            rej = 1;
        }
        fprintf(fout, "%15s    %7.3e   %d\n", loc_vec[i].c_str(),  lfdr_vec[i], rej);
    }

    fclose(fout);


}





double BAGSE_logistic::compute_log10_BF(double beta, double se_beta, double phi){

    if(se_beta == 0)
        return 0;
    int size = 4;
    double z2 = pow((beta/se_beta), 2.0);
    double v2 = pow(se_beta, 2.0);
    double w2 = pow(phi,2);
    double val = 0.5*log(v2/(v2+w2)) + 0.5*z2*(w2/(v2+w2));
    return    val/log(10);
}





BAGSE_logistic::~BAGSE_logistic(){

    if(Xd != 0)
        gsl_matrix_int_free(Xd);
    if(Xc !=0)
        gsl_matrix_free(Xc);
    if(dlevel !=0)
        gsl_vector_int_free(dlevel);

}
