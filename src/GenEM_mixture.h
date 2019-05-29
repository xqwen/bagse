#include <vector>



class GenEM_mixture{

    private:

        int K; // number of parameters
        int N; // number of units

        vector<double> wts_matrix; //  (Kx1)
        vector<vector<double> > log10_BF_matrix;  //  (NxK)
        vector<double> BF_avg;     //  (Nx1)
        vector<vector<double> > P_matrix;   //  (NxK) 

        void E_step(); // compute BF_avg and P_matrix
        void M_step(); // compute P_matrix and update wts_matrix

        void EM_init(vector<vector<double> >& BF_in, vector<double> & init_wts);
        double EM_update();    



    public:
        
       double EM_run(vector<vector<double> >& BF_in, vector<double> & init_wts, double thresh);
        vector<double> get_estimate(){
            return wts_matrix;
        }
        double compute_loglik(vector<double> & input_wts);


};


double log10_weighted_sum(vector<double> &vec, vector<double> &wts);
