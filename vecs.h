#ifndef VECTORS_OPERATIONS_MODULE_H
#define VECTORS_OPERATIONS_MODULE_H

#include <math.h>
#include <stdlib.h>

typedef double (*Func)(double);



double *get_vector(Func f, int size, double origin, double step);

//get weight of node for quadrature integration
inline double weight(int i, int n, double step)
{
    return i == 0 || i == n - 1 ? step / 2 : step;
}

double get_dot(const double *C, const double *W, int n_count, double step);

double get_norm(const double *f, int n_count, double step);

void mul(double *v, double fact, int size);

double get_diff(const double *Cn, const double *Cn_1, int n_count);

double get_relative_error(const double *Cn, Func sol, int n_count,
    double origin, double step);

#endif
