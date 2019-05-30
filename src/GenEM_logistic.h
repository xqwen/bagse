#include "GenEM.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


class GenEM_logistic : public GenEM {

    private:

        int K; // number of weight parameters
        int kc;
        int kd;
        int p; // number of logistic covariates including intercept
        int grid_size;


        vector<double> wts_matrix; //  (Kx1)
        gsl_vector *alpha_vec;  // px1

        gsl_vector *pip_vec;
        gsl_vector *prior_vec;
        vector<double> BF_avg;     //  (Nx1)
        vector<vector<double> > P_matrix;   //  Nx(K+1) 

        
        gsl_matrix *Xc; //continuous covariates 
        gsl_matrix_int *Xd; //discrete/categorical covariates    
        gsl_vector_int *dlevel;


        void EM_init(vector<vector<double> >& BF_in, gsl_vector *init_alpha, vector<double> & init_wts, gsl_matrix_int *Xd_in, gsl_vector_int *dlevel_in, gsl_matrix *Xc_in);
        double EM_update();    

    public:

        double EM_run(vector<vector<double> >& BF_in, gsl_vector *init_alpha, vector<double> & init_wts, gsl_matrix_int *Xd_in, gsl_vector_int *dlevel, gsl_matrix *Xc_in, double thresh);

        vector<double> get_wts_estimate(){
            return wts_matrix;
        }

        gsl_vector* get_alpha_estimate(){
            return  alpha_vec;
        }
        
        gsl_vector* get_pip_vec(){
            return pip_vec;
        }

        double compute_loglik(vector<double> & input_wts);
        double compute_loglik(gsl_vector * input_alpha);

        ~GenEM_logistic();


};


