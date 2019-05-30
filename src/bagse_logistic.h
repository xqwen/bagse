using namespace std;
#include <vector>
#include <map>
#include "GenEM_logistic.h"

class BAGSE_logistic {

    private:

        int N;  // unit number

        int grid_size;

        int kc; // number of continuous annotations
        int kd; // number of discrete/categorical annotations
        
        vector<vector<double> > log10_BF_matrix;
        vector<string> loc_vec;

        GenEM_logistic gem;
        double final_loglik;
        

        gsl_matrix *Xc; //continuous covariates
        gsl_matrix_int *Xd; //discrete/categorical covariates
        gsl_vector_int *dlevel;

    
        vector<string> cvar_name_vec;
        vector<string> dvar_name_vec;

        map<string, int> gene_hash;

        vector<double>  make_grid(double min, double max);
        double compute_log10_BF(double beta, double se, double phi);
        int count_factor_level(int col);

    public:

        void load_data(char *filename, int use_zval, char *annot_file);
        void load_annotation(char *annot_file);
        void run(double thresh);
        void fdr_control(char *fdr_out, double fdr_level);


        ~BAGSE_logistic();

};
