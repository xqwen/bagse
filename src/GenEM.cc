using namespace std;

#include "GenEM.h"
#include <math.h>

double GenEM::log10_weighted_sum(vector<double> &vec, vector<double> &wts){


    double max = vec[0];
    for(size_t i=0;i<vec.size();i++){
        if(vec[i]>max)
            max = vec[i];
    }
    double sum = 0;
    for(size_t i=0;i<vec.size();i++){
        sum += wts[i]*pow(10, (vec[i]-max));
    }

    return (max+log10(sum));
}


