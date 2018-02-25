#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_statistics_double.h>

/*
   Tom Wallace <twalla11@masonlive.gmu.edu>
   STAT 778
   Spring 2018
   HW 2, Problem 2

   Compilation requires GNU Scientific Library (gsl)
   
   Compilation:
       gcc -Wall -I/usr/include -c hw2_2.c
       gcc -L/usr/lib/i386-linux-gnu -o hw2_2 -lgsl -lgslcblas -lm hw2_2.o
 
   You may need to change the arguments passed to -I and -L depending on where
   the gsl header files and libgsl, respectively, live on your system

   Usage:
       Output is printed to stdout - you probably want to pipe to a text file
       Example: ./hw2_2 > output.csv
*/

void run_simulation(int n, double mu, double ss, int seed);

int main(){

	int n_iter = 1000;
	int i;

	printf("n, mean, var, SE_mean, SE_var, mean_lower95, mean_upper95, var_lower95, var_upper95\n");

	// simulation with n=50
	for(i=0; i < n_iter; i++){
		run_simulation(50, -0.5, 2.0, rand());
	}

	// simulation with n=100
	for(i=0; i < n_iter; i++){
		run_simulation(100, -0.5, 2.0, rand());
	}

	// simulation with n=200
	for(i=0; i < n_iter; i++){
		run_simulation(200, -0.5, 2.0, rand());
	}
}

void run_simulation(int n, double mu, double ss, int seed){

	// set up rng
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	// array in which to store normal RVs
	double *A = malloc(n * sizeof(double));

	// generate normal RVs and store in array
	int i;
	for(i=0; i < n; i++){
		double x = gsl_ran_gaussian(r, sqrt(ss));
		x += mu;
		A[i] = x;
	}

	// point estimate of mu
	double mean = gsl_stats_mean(A, 1, n);

	// point estimate of ss
	double var = gsl_stats_variance(A, 1, n);

	// get rid of A
	free(A);

	// SE of mean
	double se_mean = sqrt(var / n);

	// SE of ss
	double se_var = sqrt((2*pow(var, 2)) / (n-1));

	// 95% CI of mean
	double mean_u = mean + 1.96*se_mean;
	double mean_l = mean - 1.96*se_mean;

	// 95% CI of var
	double var_u = var + 1.96*se_var;
	double var_l = var - 1.96*se_var;

	printf("%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", n, mean, var, se_mean, se_var, mean_l, mean_u, var_l, var_u);
}
