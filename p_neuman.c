#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#ifndef PATH
#define PATH "graph.plt"        //file to storing data
#endif

#ifndef BAR_WIDTH
#define BAR_WIDTH 70            //width of progress bar in chars
#endif




typedef double (*Func)(double);     //type of function

int n_count = 5001;
int iters = 1000;
const double R = 20;            //limit of integration
const double b = 1;             //birth coeff
const double s = 1;             //death coeff
const double A = 20;            //kernel parameters
const double B = 0;             
double step;                    //step of nodes
#ifdef SHOUT
int sharps = 0;                 //for progress bar
#endif



//death kernel
static inline double w(double x)
{
    double ex = exp(-fabs(x));
    double Axx = A * x * x;
    
    return ex *
        (Axx/3 - 16.0/9*A*fabs(x) + 56.0/27*A + B/3) /
        (1 + ex * (Axx + B));
}



//birth kernel
static inline double m(double x)
{
    return exp(-2 * fabs(x));
}



//origin function
static inline double origin(double x)
{
    return 0.0;
}



//accurate solution
static inline double sol(double x)
{
    return exp(-fabs(x)) * (A*x*x + B);
}



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
 


//get weight of node for quadrature integration
static inline double weight(int i)
{
    return i == 0 || i == n_count - 1 ? step/2 : step;
}



//get dot product of two functions
double get_dot(const double *C, const double *W)
{
    double res = 0;
    int i;

    if(W == NULL){
        for(i = 0; i < n_count; i++){
            res += C[i]*weight(i);
        }
    }else{
        for(i = 0; i < n_count; i++){
            res += W[i] * C[i] * weight(i);
        }
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



//multiply two complex vectors storing result in the first one
static inline void comp_mul(fftw_complex *f, const fftw_complex *g, int n)
{
    double re;
    double im;
    int i;

    for(i = 0; i < n; i++){
        re = f[i][0] * g[i][0] - f[i][1] * g[i][1];
        im = f[i][0] * g[i][1] + f[i][1] * g[i][0];
        
        f[i][0] = re;
        f[i][1] = im;
    }
}



//fourier convolution for vectors F and G with the next features:
//  1. size of F=2*size, size of G=size
//  2. we pade G by zeros to extend it to size 2*size
//  3. result size is 2*size, but only points [size, 2*size] matters for us
double *fourier_conv(double *F, double *G, int size)
{
    double *g1 = malloc(sizeof(double) * 2 * size);
    fftw_complex *tmp1 = fftw_alloc_complex(size + 1);
    fftw_complex *tmp2 = fftw_alloc_complex(size + 1);
    double *result = malloc(sizeof(double) * 2 * size);
    int i;

    for(i = 0; i < size; i++){
        g1[i] = G[i];
    }

    for(; i < 2*size; i++){
        g1[i] = 0.0;
    }

    fftw_plan pF;
    fftw_plan pf2;
    fftw_plan pb;

    pF = fftw_plan_dft_r2c_1d(size * 2, F, tmp1, FFTW_ESTIMATE);
    pf2 = fftw_plan_dft_r2c_1d(size * 2, g1, tmp2, FFTW_ESTIMATE);
    pb = fftw_plan_dft_c2r_1d(size * 2, tmp1, result, FFTW_ESTIMATE);

    fftw_execute(pF);
    fftw_execute(pf2);
    comp_mul(tmp1, tmp2, size + 1);
    fftw_execute(pb);

    fftw_destroy_plan(pF);
    fftw_destroy_plan(pf2);
    fftw_destroy_plan(pb);
    
    fftw_free(tmp2);
    fftw_free(tmp1);
    free(g1);
    mul(result, step/(size * 2), size * 2);

    return result;
}



//get C-norm of difference of two function
double get_diff(const double *Cn, const double *Cn_1)
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
double get_relative_error(const double *Cn, Func sol)
{
    double curr;
    double max = 0;
    int i;
    double x = -R;
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



//store calculated solution
void store_solution(const double *C)
{
    FILE *out = fopen(PATH, "w");
    double x = -R;
    double f;
    int i;

    for(i = 0; i < n_count; i++){
        f = sol(x);
        fprintf(out,
            "%15.7lf %15.7lf %15.7lf %15.7lf\n",
            x,
            C[i] + 1,
            f + 1,
            fabs(f - C[i]) / (f + 1)
        );
        x += step;
    }

    fclose(out);
}



//iterate process
//Nota bene: M size is 2*n_count (it's extended)
void twin_iterate(double **Cn, double **Cn_1, double *M, double *W,
    double N, double nw)
{
    double *conv;
    double bswx;
    int i;

    conv = *Cn;
    *Cn = *Cn_1;
    *Cn_1 = conv;

    conv = fourier_conv(M, *Cn_1, n_count);

    for(i = 0; i < n_count; i++){
        bswx = b + s*W[i];
        (*Cn)[i] = (b*conv[i + n_count - 1] + M[i+n_count/2]*N +
            s*(M[i+n_count/2]*nw - W[i])) / bswx;
    }

    free(conv);
}



//get solution of twin_equation of current parameter N
double *get_solution(double N, double *M, double *W, double nw)
{
    double *Cn = malloc(sizeof(double) * n_count);  //current iteration
    double *Cn_1 = malloc(sizeof(double)*n_count);  //the last iteration
    int i;

    Cn = get_vector(&origin, n_count, -R, step);;

#   ifdef SHOUT
    printf("Progress: ");
    fflush(stdout);
#   endif

    for(i = 0; i < iters; i++){
        //computing
        twin_iterate(&Cn, &Cn_1, M, W, N, nw);

        //working with interface
#       ifndef SHOUT
        printf("Iteration: %i\nDifference: %40.30f\n", i,
            get_diff(Cn, Cn_1));
#       else
        while(sharps < (double)BAR_WIDTH * i / iters){
            putchar('#');
            sharps++;
        }
        fflush(stdout);
#       endif
    }

#   ifdef SHOUT
    for(;sharps < BAR_WIDTH; sharps++){
        putchar('#');
    }
    putchar('\n');
    sharps = 0;
#   endif

    free(Cn_1);

    return Cn;
}



int str2int(const char *str)
{
    int i;
    int res = 0;
    
    for(i = 0; str[i]; i++){
        res *= 10;
        res += str[i] - '0';
    }

    return res;
}



void set_values(int argc, char **argv)
{
    if(argc < 2){
        return;
    }

    iters = str2int(argv[1]);

    if(argc < 3){
        return;
    }

    n_count = str2int(argv[2]);
}



int main(int argc, char **argv)
{
    set_values(argc, argv);
    step = 2 * R / (n_count - 1);

    double *W = get_vector(&w, n_count, -R, step);      //kernel w
    double *M = get_vector(&m, 2*n_count, -2*R, step);  //kernel m
    double *solution;
    double *h;
    double *g;
    double N;
    double nw = get_dot(W, NULL);

    printf("***Getting g(x)...\n");
    g = get_solution(0, M, W, nw);

    printf("***Getting h(x)...\n");
    h = get_solution(1, M, W, nw);

    printf("***Getting solution...\n");
    N = s * get_dot(W, g);
    N = N / (1 - s * get_dot(W, h) + N);
    solution = get_solution(N, M, W, nw);

    printf("\nError: %60.40lf\n", get_relative_error(solution, &sol));
    store_solution(solution);

    free(M);
    free(W);
    free(h);
    free(g);
    free(solution);

    return 0;
}
