using namespace std;

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "bagse_mixture.h"

#define NINF -999999

void BAGSE_mixture::load_data(char *filename, int use_zval){



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


    category_map[string("0")] = 0;
    category_rmap[0] = string("0");
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
        
        if(!(ins>>annot)){
            fprintf(stderr, "\nError: unexpected data format in input\n\n");
            exit(1);
        };
        loc_vec.push_back(loc_id);
        beta_vec.push_back(beta);
        se_vec.push_back(se_beta);

        if(se_beta<min)
            min = se_beta;

        if(beta*beta - se_beta*se_beta > max){
            max = pow(beta,2)-pow(se_beta,2);
        }



        if(category_map.find(annot) == category_map.end()){
            int cat = category_map.size();

            category_map[annot] = cat;
            category_rmap[cat] = annot;
        }
        annot_vec.push_back(category_map[annot]);

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


    annot_size = category_map.size();

    K = annot_size*(grid_size+1);

    N = beta_vec.size();


    fprintf(stderr, "Initializing ... \n");
    //fprintf(stderr, "N=%d\t K=%d\t L=%d\t  M=%d\n", N, grid_size, annot_size, K);

    // compute BF factors
    for(int i=0;i<N;i++){
        vector<double> bf_vec(K,NINF);
        for(int j = 0; j< grid_size;j++){
            
            int index = (1+grid_size)*annot_vec[i]+j;
            double log10_bf = 0;
            if(j>0){
                log10_bf = compute_log10_BF(beta_vec[i], se_vec[i],grid_vec[j-1]);
            }
            //printf("index = %d log10_bf = %f %d %d\n", index, log10_bf, annot_vec[i], j);
            bf_vec[index] = log10_bf;
        }


       log10_BF_matrix.push_back(bf_vec);  
    }

    /*
    for(int i=0;i<N;i++){
        for(int j=0; j<K;j++){
            printf("%7.3f ", log10_BF_matrix[i][j]);
        }
        printf("\n");
    }
    */

}


vector<double> BAGSE_mixture::make_grid(double phi_min, double phi_max){

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




void BAGSE_mixture::run(double thresh){  

    // initial setting
    vector<double> wts_vec(K, (1.0/N)/(grid_size*annot_size) );
    for(int i=0;i<K;i+= grid_size+1){
        wts_vec[i] = (1-1.0/N)/(annot_size);
    }
    
    /*
    double sum = 0;
    for(int i=0;i<K;i++){
        sum += wts_vec[i];
    }

    printf("sum = %f\n",sum);
    exit(0);
    */
    fprintf(stderr, "Running EM ... \n\n");
    final_loglik = gem.EM_run(log10_BF_matrix, wts_vec, thresh);
    vector<double> est_vec = gem.get_estimate();
    vector<double> alpha_vec;

    for(int cat=0;cat<category_map.size();cat++){
        double prob1 = 0;
        for(int j=0;j<grid_size;j++){
            int index = cat*(grid_size+1)+j+1;
            prob1+= est_vec[index];
        }
        double prob0 = est_vec[cat*(grid_size+1)];
        alpha_vec.push_back(log(prob1/prob0));
    }
    fprintf(stderr, "\n\n");


    vector<double> ci_vec = find_CI(est_vec, 0,alpha_vec[0],  0);
    printf("%15s   %10s   %7.3f    %7.3f %7.3f\n", "Intercept", "0", alpha_vec[0], ci_vec[0], ci_vec[1]);

    for(int cat =1; cat <category_map.size(); cat++){
        ci_vec = find_CI(est_vec, cat, alpha_vec[cat]-alpha_vec[0],  alpha_vec[0]);
        char header[128];
        sprintf(header, "annot.%d", cat);
        printf("%15s   %10s   %7.3f    %7.3f %7.3f\n",header, category_rmap[cat].c_str(),alpha_vec[cat]-alpha_vec[0], ci_vec[0], ci_vec[1]);
    }
    
    fprintf(stderr, "\n\n");

}



vector<double> BAGSE_mixture::find_CI(vector<double> & est_vec, int cat, double value, double contrast){

    double total_wts = 0;

    double pi0 = est_vec[(grid_size+1)*cat];
    double pi1 = 0;
    for(int j=0;j<grid_size;j++){
        int index = cat*(grid_size+1)+j+1;
        pi1+= est_vec[index];
    }   


    vector<double> input_vec = est_vec;
    input_vec[(grid_size+1)*cat] = (pi0+pi1)/(1+exp(contrast));
    double npi1 = (pi0+pi1)*exp(contrast)/(1+exp(contrast));
    for(int j=0;j<grid_size;j++){
        int index = cat*(grid_size+1)+j+1;
        input_vec[index] = npi1*est_vec[index]/pi1;
    }

    double loglik = gem.compute_loglik(input_vec);
    double se = sqrt(pow(value,2)/(2*(final_loglik - loglik)));
    vector<double> rstv;
    rstv.push_back(value - 1.96*se);
    rstv.push_back(value + 1.96*se);
    return rstv;

}


void BAGSE_mixture::fdr_control(char *fdr_out, double fdr_level){

    vector<vector<double> > P_matrix = gem.get_P_matrix();
    vector<double> lfdr_vec;
    for(int i=0;i<N;i++){
        int index = annot_vec[i]*(grid_size+1);
        lfdr_vec.push_back(P_matrix[i][index]);
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
        fprintf(fout, "%15s  %10s   %7.3e   %d\n", loc_vec[i].c_str(), category_rmap[annot_vec[i]].c_str(),  lfdr_vec[i], rej);
    }

    fclose(fout);


}





double BAGSE_mixture::compute_log10_BF(double beta, double se_beta, double phi){

    if(se_beta == 0)
        return 0;
    int size = 4;
    double z2 = pow((beta/se_beta), 2.0);
    double v2 = pow(se_beta, 2.0);
    double w2 = pow(phi,2);
    double val = 0.5*log(v2/(v2+w2)) + 0.5*z2*(w2/(v2+w2));
    return    val/log(10);
}

