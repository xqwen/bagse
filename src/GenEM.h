#ifndef _GENEM_H_
#define _GENEM_H_

#include <vector>

class GenEM {

    protected:
        int N;  
        vector<vector<double> > log10_BF_matrix;
        
        virtual double EM_update() = 0;
        double log10_weighted_sum(vector<double> &vec, vector<double> &wts);
};
#endif
