using namespace std;
#include <vector>
#include <map>
#include "GenEM_mixture.h"

class BAGSE {

    private:
        int N;  // unit number
        int K; // grid number
         
        int grid_size;
        int annot_size;

        
        map<string, int> category_map;
        map<int, string> category_rmap;
        vector<int> annot_vec;
        vector<vector<double> > log10_BF_matrix;
        vector<string> loc_vec;

        GenEM_mixture gem;
        double final_loglik;
        

        vector<double>  make_grid(double min, double max);
        double compute_log10_BF(double beta, double se, double phi);
        vector<double> find_CI(vector<double> &vec, int cat, double value, double contrast);

    public:
        void run(double thresh);
        void load_data(char *filename, int use_zval);
        void fdr_control(char *fdr_out, double fdr_level);
};
