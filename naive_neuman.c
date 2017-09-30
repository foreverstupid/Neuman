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



//get convolution of two functions
double *convolve(double *f, Func g)
{
	double *res = malloc(sizeof(double) * N_COUNT);
	int i;
	int j;
	double x = -R;
	double y;

	for(i = 0; i < N_COUNT; i++){
		y = -R;
		res[i] = 0.0;
		for(j = 0; j < N_COUNT; j++){
			res[i] += weight(j) * f[j] * g(x -y); 
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

	conv = convolve(*Cn_1, &m);

	for(i = 0; i < N_COUNT; i++){
		bswx = b + s*W[i];
		(*Cn)[i] = (b*conv[i] + M[i]*N +
			s*(M[i]*nw - W[i])) / bswx;
	}

	free(conv);
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
	double *h;
	double *g;
	double N;
	double nw = get_dot(W, NULL);

	printf("***Getting g(x)...\n");
	g = get_solution(0, M, W, nw);

	printf("***Getting h(x)...\n");
	h = get_solution(1, M, W, nw);

	printf("***Getting solution...\n");
	N = get_dot(W, g);
	N = -N / (get_dot(W, h) - 1 - N);
	solution = get_solution(N, M, W, nw);

	printf("\nError: %20.10lf\n", get_relative_error(solution, &sol));
	store_solution(solution);

	free(M);
	free(W);
	free(h);
	free(g);
	free(solution);

	return 0;
}
