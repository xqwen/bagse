using namespace std;

#include "GenEM_mixture.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void GenEM_mixture::EM_init(vector<vector<double> >& log10_BF_in, vector<double> & init_wts){

    K = int(init_wts.size());
    N = int(log10_BF_in.size());


    wts_matrix = init_wts;
    log10_BF_matrix = log10_BF_in;
    P_matrix = log10_BF_in; 
 
    /*
    double sum = 0;
    for(int i=0;i<K;i++){
        printf("%7.3e %7.3f", wts_matrix[i], log10(wts_matrix[i]));
        sum += wts_matrix[i];
    }
    printf("\n%7.3f\n",sum);
    */

}





double GenEM_mixture::EM_update(){
    // compute BF_avg
    //
    // E-step
    double log10_lik = 0; 
    vector<double> log10_BF_avg(N,0);

    for(int i=0;i<N;i++){
        log10_BF_avg[i] = log10_weighted_sum(log10_BF_matrix[i], wts_matrix);
        log10_lik += log10_BF_avg[i];
    }

    vector<double> new_wts(K,0);
    for(int i=0;i<N;i++){
        for(int j=0;j<K;j++){
            double val = log10_BF_matrix[i][j] + log10(wts_matrix[j]) - log10_BF_avg[i];
            P_matrix[i][j] = pow(10, val);
            new_wts[j] += P_matrix[i][j];
        }
    }

    for(int j=0;j<K;j++){
        wts_matrix[j] = new_wts[j]/N;
    }

    return log10_lik;
}


double GenEM_mixture::compute_loglik(vector<double> &input_wts){

    double log10_lik = 0;
    for(int i=0;i<N;i++){
        log10_lik += log10_weighted_sum(log10_BF_matrix[i], input_wts);
    }
    return log10_lik/log10(exp(1));
}



double GenEM_mixture::EM_run(vector<vector<double> >& log10_BF_in, vector<double> & init_wts, double thresh){

    EM_init(log10_BF_in, init_wts);
    double log10_lik = EM_update();
    int iter = 1;
    while(1){
        fprintf(stderr,"Iter %d\t\t%7.3f\n", iter, log10_lik/log10(exp(1)));
        double log10_lik_new = EM_update();
        if(log10_lik_new - log10_lik <= log10(exp(1))*thresh){
            break;
        }else{
            log10_lik = log10_lik_new;
            iter++;
        }
    }
    return log10_lik/log10(exp(1));
}



