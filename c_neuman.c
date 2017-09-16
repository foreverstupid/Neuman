#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/time.h>

typedef double (*func)(double);

#ifndef N_COUNT
#define N_COUNT 5001 	//count of nodes
#endif

#ifndef ITERS
#define ITERS 1000 		//iterations
#endif

#ifndef PATH
#define PATH "graph.plt"	//path to data file
#endif

const double p = 1;		//equation parameters
const double b = 1;
const double s = 1;
#ifdef EXP
const double A = 2;
const double B = 1;
#else
const double A = 1;
#endif
double N; 
const double R = 10.0;

double step;				//solution step


//birth kernel
static inline double m(double x)
{
#	ifndef EXP
	return p / (M_PI * (x * x + p * p));
#	else
	return exp(-2 * fabs(x));
#	endif
}



//death kernel
static inline double w(double x)
{
#	ifndef EXP
	return A / (x * x + 9 * p * p);
#	else
	double ex = exp(-fabs(x));
	double Axx = A * x * x;
	
	return ex *
		(Axx/3 - 16.0/9*A*fabs(x) + 56.0/27*A + B/3) /
		(1 + ex * (Axx + B));
#	endif
}



//function f from equation g = Kg + f
static inline double f(double x)
{
	double swx = s * w(x);

	return (N * m(x) - swx) / (b + swx);
}



//integral operator kernel
static inline double kernel(double x, double y)
{
	return b * m(x - y) / (b + s * w(x));
}



//accurate solution
static inline double solution(double x)
//Nota bene: we replace variable function: C(x) = solution(x) + 1
{
#	ifndef EXP
	double xx = x * x;
	return 24.0 / 71.0 / (xx + 1) + 40.0 / 71.0 / (xx + 4);
#	else
	return exp(-fabs(x)) * (A*x*x + B);
#	endif
}



//create vector of function values in nodes
//Nota bene: we use parity of functions
double *make_vector(func f)
{
	int i;
	double *res = malloc(sizeof(double) * N_COUNT);
	double x = -R;

	for(i = 0; i < N_COUNT / 2 + 1; i++){
		res[i] = res[N_COUNT - i - 1] = f(x);
		x += step;
	}

	return res;
}



//weights for quadrature formula
static inline double weight(int j)
{
	return j == 0 || j == N_COUNT - 1  ? step / 2 : step;
}



//create matrix from integral operator
//Nota bene: we hold matrix in one-dimensional array
double *make_matrix()
{
	int i;
	int j;
	double *K_mat = malloc(sizeof(double) * N_COUNT * N_COUNT);

	double x = -R;
	double y;

	for(i = 0; i < N_COUNT; i++){
		y = -R;
		for(j = 0; j < N_COUNT; j++){
			K_mat[i * N_COUNT + j] = weight(j) * kernel(x, y);
			y += step;
		}
		x += step;
	}

	return K_mat;
}



//get relative error of current solution
double get_err(double *C)
//Nota bene: error is calculated consider that C(x) = solution(x) + 1
//           and that C(x) = C(-x)
{
	int i;
	double res = 0;
	double err;
	double x = -R;
	double sx;

	for(i = 0; i < N_COUNT / 2 + 1; i++){
		sx = solution(x);
		err = fabs(C[i] - sx) / (sx + 1);
		if(err > res){
			res = err;
		}

		x += step;
	}

	return res;
}



//calculating current iteration
//C_(n+1) = K * C_n + f
void calculate(double **C, double *F, double *K)
{
	int i;
	int j;
	int i_shr;
	double x;
	double *res = malloc(sizeof(double) * N_COUNT);

	for(i = 0; i < N_COUNT; i++){
		x = F[i];
		i_shr = i * N_COUNT;
		for(j = 0; j < N_COUNT; j++){
			x += K[i_shr + j] * (*C)[j];
		}
		res[i] = x;
	}

	free(*C);
	*C = res;
}



//save data to data file
void push_data(double *C)
//Nota bene: remember that C(x) = solution(x) + 1
{
	int i;
	FILE *out = fopen(PATH, "w");
	double x = -R;
	double sx;
	
	for(i = 0; i < N_COUNT; i++){
		sx = solution(x);
		fprintf(
			out,
			"%12.5lf%12.5lf%12.5lf%12.5lf\n",
			x,
			C[i] + 1,
			sx + 1,
			fabs(C[i] - sx) / (sx + 1)
		);
		x += step;
	}

	fclose(out);
}



//get time of calculating
//(we don't think about preparation for example matrix constructing)
long get_time(struct timeval *tm1, struct timeval *tm2)
{
	long microseconds = tm2->tv_usec - tm1->tv_usec;
	long seconds = tm2->tv_sec - tm1->tv_sec;

	if(microseconds < 0){
		seconds--;
		microseconds += 1000000;
	}

	return seconds * 1000 + microseconds / 1000;
}




int main()
{
	struct timeval tm1, tm2;
#	ifndef EXP
	N = A*M_PI/p * (A+5*p*p) * (A+8*p*p) / (A*A + 21*A*p*p + 120*p*p);
#	else
	N = 2.0/3*B + 52.0/27*A;
#	endif
	step = 2.0 * R / (N_COUNT - 1);
	
	int iters = 0;
	double *F = make_vector(&f);
	double *C = make_vector(&f);
	double *K = make_matrix();

	gettimeofday(&tm1, NULL);
	while(iters < ITERS){
		calculate(&C, F, K);
		iters++;

#		ifndef SHOUT
		printf(
			"Iteration: %i\n"
			"Error: %20.15lf\n",
			iters,
			get_err(C)
		);
#		endif
	}

	gettimeofday(&tm2, NULL);
	printf("Time: %li milliseconds\n", get_time(&tm1, &tm2));
	printf("Error: %12.7lf\n", get_err(C));

	push_data(C);

	free(F);
	free(C);
	free(K);

	return 0;
}
