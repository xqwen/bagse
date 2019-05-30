using namespace std;

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "GenEM_logistic.h"
#include "logistic.h"

void GenEM_logistic::EM_init(vector<vector<double> >& log10_BF_in, gsl_vector *init_alpha, vector<double> & init_wts, gsl_matrix_int *Xd_in, gsl_vector_int *dlevel_in, gsl_matrix *Xc_in){


    Xc = 0;
    Xd = 0;
    dlevel = 0;
    alpha_vec = 0;
    pip_vec = 0;
    prior_vec = 0;

    K = int(init_wts.size());
    N = int(log10_BF_in.size());
    
    if(Xc_in != 0)
        kc = Xc_in->size2;
    else 
        kc=0;

    if(Xd_in !=0)
        kd = Xd_in->size2;
    else
        kd = 0;

    p = init_alpha->size;

    log10_BF_matrix = log10_BF_in;
    P_matrix = log10_BF_in; 


    alpha_vec = gsl_vector_calloc(p);
    gsl_vector_memcpy(alpha_vec, init_alpha);
    
    
    wts_matrix = init_wts;
    grid_size = wts_matrix.size();
    double sum = 0;
    for(int i=0;i<grid_size;i++){
        sum += wts_matrix[i];
    }



    if(kc>0){
        Xc = gsl_matrix_calloc(N,kc);
        gsl_matrix_memcpy(Xc, Xc_in);
    }
    if(kd>0){
        Xd = gsl_matrix_int_calloc(N,kd);
        gsl_matrix_int_memcpy(Xd, Xd_in);
        dlevel = gsl_vector_int_calloc(kd);
        gsl_vector_int_memcpy(dlevel, dlevel_in);
    }

    pip_vec = gsl_vector_calloc(N);
    prior_vec = gsl_vector_calloc(N);
    // init
    for(int i=0;i<N;i++){
        gsl_vector_set(prior_vec,i,1e-3);
    }   

    
}


double GenEM_logistic::compute_loglik(gsl_vector *alpha_input){
     
    gsl_vector * new_prior_vec = gsl_vector_calloc(N);
     if(kc==0 && kd!=0){
        logistic_cat_pred(alpha_input, Xd, dlevel, new_prior_vec);
    }


    if(kc!=0 && kd==0){
        logistic_cont_pred(alpha_input,Xc,new_prior_vec);
    }

    if(kc!=0 && kd!=0){
        logistic_mixed_pred(alpha_input, Xd, dlevel, Xc, new_prior_vec);
    }


    double log10_lik = 0;
    vector<double> log10_BF_avg(N,0);
    for(int i=0;i<N;i++){

        vector<double> gene_wts;
        double pi1 = gsl_vector_get(new_prior_vec,i);
        gene_wts.push_back(1-pi1);
        for(int j=0;j<grid_size;j++){
            gene_wts.push_back(pi1*wts_matrix[j]);
        }
        log10_BF_avg[i] = log10_weighted_sum(log10_BF_matrix[i], gene_wts);
        log10_lik += log10_BF_avg[i];
    }

    gsl_vector_free(new_prior_vec);

    return log10_lik/log10(exp(1));

}




double GenEM_logistic::EM_update(){

    double log10_lik = 0; 
    vector<double> log10_BF_avg(N,0);
    
    vector<vector<double> > prior_matrix;
    for(int i=0;i<N;i++){
        
        vector<double> gene_wts;
        double pi1 = gsl_vector_get(prior_vec,i);
        gene_wts.push_back(1-pi1);
        for(int j=0;j<grid_size;j++){
            gene_wts.push_back(pi1*wts_matrix[j]);
        }
        log10_BF_avg[i] = log10_weighted_sum(log10_BF_matrix[i], gene_wts);
        log10_lik += log10_BF_avg[i];
        prior_matrix.push_back(gene_wts);
    }    


     vector<double> new_wts(grid_size,0);
     double psum = 0;
     for(int i=0;i<N;i++){
        for(int j=0;j<grid_size+1;j++){
            double val = log10_BF_matrix[i][j] + log10(prior_matrix[i][j]) - log10_BF_avg[i];
            P_matrix[i][j] = pow(10, val);
            if(j>0){
                 new_wts[j-1] += P_matrix[i][j];
                 psum += P_matrix[i][j];
            }
        } 
        gsl_vector_set(pip_vec,i,1-P_matrix[i][0]); 
     }
        
        
     // update alpha_vec and pip_vec
    if(kc==0 && kd!=0){

        logistic_cat_fit(alpha_vec, Xd, dlevel,pip_vec,0,0);
        logistic_cat_pred(alpha_vec, Xd, dlevel, prior_vec);
    }


    if(kc!=0 && kd==0){
        logistic_cont_fit(alpha_vec, Xc, pip_vec,0,0);
        logistic_cont_pred(alpha_vec,Xc,prior_vec);
    }

    if(kc!=0 && kd!=0){
        logistic_mixed_fit(alpha_vec, Xd, dlevel, Xc, pip_vec, 0,0);
        logistic_mixed_pred(alpha_vec, Xd, dlevel, Xc, prior_vec);
    }




    for(int j=0;j<grid_size;j++){
        wts_matrix[j] = new_wts[j]/psum;
    }
    return log10_lik;



}






double GenEM_logistic::EM_run(vector<vector<double> >& log10_BF_in, gsl_vector *init_alpha, vector<double> & init_wts, gsl_matrix_int *Xd_in, gsl_vector_int *dlevel_in, gsl_matrix *Xc_in, double thresh){

    EM_init(log10_BF_in, init_alpha, init_wts, Xd_in, dlevel_in, Xc_in);
    double log10_lik = EM_update();
    int iter = 1;
    while(1){
        fprintf(stderr,"Iter %d\t\t%7.3f ", iter, log10_lik/log10(exp(1)));
        /*       
                 for(int i=0;i<alpha_vec->size;i++){
                 fprintf(stderr, "%7.3f  ", gsl_vector_get(alpha_vec,i));
                 }
                 */
        fprintf(stderr, "\n");

        double log10_lik_new = EM_update();
        if(log10_lik_new - log10_lik <= log10(exp(1))*thresh){
            break;
        }else{
            log10_lik = log10_lik_new;
            iter++;
        }
    }
    fprintf(stderr, "\n\n");
    return log10_lik/log10(exp(1));
}




GenEM_logistic::~GenEM_logistic(){

    if(Xd != 0)
        gsl_matrix_int_free(Xd);
    if(Xc !=0)
        gsl_matrix_free(Xc);
    if(dlevel !=0)
        gsl_vector_int_free(dlevel);

    if(pip_vec!=0)
        gsl_vector_free(pip_vec);

    if(prior_vec!=0)
        gsl_vector_free(prior_vec);

    if(alpha_vec!=0)
        gsl_vector_free(alpha_vec);



}





