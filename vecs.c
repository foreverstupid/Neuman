#include "vecs.h"

//make vector from scalar function
double *get_vector(Func f, int size, double origin, double step)
{
    double *res = malloc(sizeof(double) * size);
    double x = origin;
    int i;
    
    for(i = 0; i < size; i++){
        res[i] = f(x);
        x += step;
    }

    return res;
}



//get dot product of two functions
double get_dot(const double *C, const double *W, int n_count, double step)
{
    double res = 0;
    int i;

	for(i = 0; i < n_count; i++){
		res += W[i] * C[i] * weight(i, n_count, step);
	}

    return res;
}



//get 1-norm of function
double get_norm(const double *f, int n_count, double step)
{
    double res = 0;
    int i;

	for(i = 0; i < n_count; i++){
		res += f[i] * weight(i, n_count, step);
	}

    return res;
}



//multiply vector by factor
void mul(double *v, double fact, int size)
{
    int i;

    for(i = 0; i < size; i++){
        v[i] *= fact;
    }
}



//get C-norm of difference of two function
double get_diff(const double *Cn, const double *Cn_1, int n_count)
{
    double curr;
    double max = 0;
    int i;

    for(i = 0; i < n_count; i++){
        curr = fabs(Cn[i] - Cn_1[i]);
        if(curr > max){
            max = curr;
        }
    }

    return max;
}



//get relative C-norm of error
double get_relative_error(const double *Cn, Func sol, int n_count,
	double origin, double step)
{
    double curr;
    double max = 0;
    int i;
    double x = origin;
    double f;

    for(i = 0; i < n_count; i++){
        f = sol(x);
        curr = fabs(Cn[i] - f) / (f + 1);
        if(curr > max){
            max = curr;
        }
        x += step;
    }

    return max;
}
