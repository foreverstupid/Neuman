#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#ifndef PATH
#define PATH "graph.plt"		//file to storing data
#endif

#ifndef N_COUNT
#define N_COUNT 5001 			//node count	
#endif

#ifndef ITERS
#define ITERS 1000 				//iterations of solving
#endif

#ifndef BAR_WIDTH
#define BAR_WIDTH 70			//width of progress bar in chars
#endif




typedef double (*Func)(double);		//type of function

const double R = 20;			//limit of integration
const double b = 1;				//birth coeff
const double s = 1;				//death coeff
const double A = 2;				//kernel parameters
const double B = 1;				
double step;					//step of nodes
#ifdef SHOUT
int sharps = 0;					//for progress bar
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



//accurate solution
static inline double sol(double x)
{
	return exp(-fabs(x)) * (A*x*x + B);
}



//make vector from scalar function
double *get_vector(Func f)
{
	double *res = malloc(sizeof(double) * N_COUNT);
	double x = -R;
	int i;
	
	for(i = 0; i < N_COUNT; i++){
		res[i] = f(x);
		x += step;
	}

	return res;
}



//initialize C_0
void init_first_iteration(double *C_0)
{
	double x = -R;
	int i = 0;

	for(i = 0; i < N_COUNT; i++){
		C_0[i] = 0;
		x += step;
	}
}



//get weight of node for quadrature integration
static inline double weight(int i)
{
	return i == 0 || i == N_COUNT - 1 ? step/2 : step;
}



//get dot product of two functions
double get_dot(const double *C, const double *W)
{
	double res = 0;
	int i;

	if(W == NULL){
		for(i = 0; i < N_COUNT; i++){
			res += C[i]*weight(i);
		}
	}else{
		for(i = 0; i < N_COUNT; i++){
			res += W[i] * C[i] * weight(i);
		}
	}

	return res;
}



//multiply vector by factor
void mul(double *v, double fact)
{
	int i;

	for(i = 0; i < N_COUNT; i++){
		v[i] *= fact;
	}
}



//multiply two complex vectors storing result in the first one
static inline void comp_mul(fftw_complex *f, const fftw_complex *g)
{
	double re;
	double im;
	int i;

	for(i = 0; i < N_COUNT; i++){
		re = f[i][0] * g[i][0] - f[i][1] * g[i][1];
		im = f[i][0] * g[i][1] + f[i][1] * g[i][0];
		
		f[i][0] = re;
		f[i][1] = im;
	}
}



//get convolution of two functions
double *convolve(double *f, double *g)
{
	fftw_complex *tmp1 = fftw_alloc_complex(N_COUNT);
	fftw_complex *tmp2 = fftw_alloc_complex(N_COUNT);
	double *result = malloc(sizeof(double) * N_COUNT);
	fftw_plan pf1;
	fftw_plan pf2;
	fftw_plan pb;

	pf1 = fftw_plan_dft_r2c_1d(N_COUNT, f, tmp1, FFTW_ESTIMATE);
	pf2 = fftw_plan_dft_r2c_1d(N_COUNT, g, tmp2, FFTW_ESTIMATE);
	pb = fftw_plan_dft_c2r_1d(N_COUNT, tmp1, result, FFTW_ESTIMATE);

	fftw_execute(pf1);
	fftw_execute(pf2);
	comp_mul(tmp1, tmp2);
	fftw_execute(pb);

	fftw_destroy_plan(pf1);
	fftw_destroy_plan(pf2);
	fftw_destroy_plan(pb);
	
	fftw_free(tmp2);
	fftw_free(tmp1);
	mul(result, step/N_COUNT);

	return result;
}



/*return n*n convolution */
double *true_conv(double *f, double *M)
{
	int i;
	int j;
	double *res = malloc(sizeof(double) * N_COUNT);
	double x = -R;
	double y;

	for(i = 0; i < N_COUNT; i++){
		y = -R;
		res[i] = 0;
		for(j = 0; j < N_COUNT; j++){
			res[i] += weight(j) * f[j] * m(x - y);
			y += step;
		}
		x += step;
	}

	return res;
}



void twin_iterate(double **Cn, double **Cn_1, double *M, double *W,
	double N, double nw)
{
	double *conv;
	double bswx;
	int i;

	conv = *Cn;
	*Cn = *Cn_1;
	*Cn_1 = conv;

	conv = convolve(*Cn_1, M);

	for(i = 0; i < N_COUNT; i++){
		bswx = b + s*W[i];
		(*Cn)[i] = (b*conv[(i + N_COUNT/2)%N_COUNT] + M[i]*N -
			s*W[i]) / bswx;
	}

	free(conv);
}



//get C-norm of difference of two function
double get_diff(const double *Cn, const double *Cn_1)
{
	double curr;
	double max = 0;
	int i;

	for(i = 0; i < N_COUNT; i++){
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

	for(i = 0; i < N_COUNT; i++){
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

	for(i = 0; i < N_COUNT; i++){
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



//get solution of twin_equation of current parameter N
double *get_solution(double N, double *M, double *W, double nw)
{
	double *Cn = malloc(sizeof(double) * N_COUNT);	//current iteration
	double *Cn_1 = malloc(sizeof(double)*N_COUNT);	//the last iteration
	int i;

	init_first_iteration(Cn);

#	ifdef SHOUT
	printf("Progress: ");
	fflush(stdout);
#	endif

	for(i = 0; i < ITERS; i++){
		//computing
		twin_iterate(&Cn, &Cn_1, M, W, N, nw);

		//working with interface
#		ifndef SHOUT
		printf("Iteration: %i\nDifference: %40.30f\n", i,
			get_relative_error(Cn, &sol));
#		else
		while(sharps < (double)BAR_WIDTH * i / ITERS){
			putchar('#');
			sharps++;
		}
		fflush(stdout);
#		endif
	}

#	ifdef SHOUT
	for(;sharps < BAR_WIDTH; sharps++){
		putchar('#');
	}
	putchar('\n');
	sharps = 0;
#	endif

	free(Cn_1);

	return Cn;
}



int main()
{
	step = 2 * R / (N_COUNT - 1);
	double *W = get_vector(&w);			//kernel w in vector form
	double *M = get_vector(&m);			//kernel m in vector form
	double *solution;
	double nw = get_dot(W, NULL);

	solution = get_solution(2.0/3*B + 52.0/27*A, M, W, nw);

	printf("\nError: %20.10lf\n", get_relative_error(solution, &sol));
	store_solution(solution);

	free(M);
	free(W);
	free(solution);

	return 0;
}
