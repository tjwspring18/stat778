#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_statistics_double.h>
#include<gsl/gsl_cdf.h>

/*
Tom Wallace <twalla11@masonlive.gmu.edu>
STAT 778
Spring 2018
Midterm

Requires GNU Scientific Library (gsl)

Compilation:
gcc -I/usr/include -c midterm.c
gcc -L/usr/lib/i386-linux-gnu midterm.o -o midterm -lgsl -lgslcblas -lm 

You may need to change the arguments passed to -I and -L depending on where
the gsl header files and libgsl, respectively, live on your system

Usage:
No arguments required, prints to stdout
*/

void generate_normal_vector(double A[], int n, double mu, double var, int seed);
void generate_exp_vector(double A[], int n, double mu, int seed);
void generate_contaminated_normal_vector(double A[], int n, double mu, double var, int seed);
double t_statistic(double A[], double B[], int n);
double t_df(double A[], double B[], int n);
double t_sig(double t_statistic, double df);
struct DF{
	double *Observations; 
	double *Source;
	double *Rank;
};
void populate_df(struct DF df, double A[], double B[], int n);
void bsort(struct DF df, int n);
void assign_rank(struct DF df, int n);
double u_statistic(struct DF df, double source, int n);
double u_norm_approx(double u, int n);

int main(){

	//number in each group
	int n = 100;

	//initialize random number generator
	time_t t;
	srand((unsigned) time(&t));

	//generate and stash RVs in respective arrays
	//we will use these for t-test
	double *A = malloc(n * sizeof(double));
	double *B = malloc(n * sizeof(double));
	generate_normal_vector(A, n, 0, 1, rand());
	generate_normal_vector(B, n, 0, 1, rand());
	
	//initialize data structure
	//copy data from vectors to data structure
	//we will use this for Wilcoxon rank-sum test
	struct DF df;
	df.Observations = (double *)malloc(sizeof(double) * 2 * n);
	df.Source = (double *)malloc(sizeof(double) * 2 * n);
	df.Rank = (double *)malloc(sizeof(double) * 2 * n);
	populate_df(df, A, B, n);

	//sort Observations and Source in ascending order
	bsort(df, n);

	//assign Rank
	assign_rank(df, n);
	
	//calculate u statistic
	//should be min(u_A, u_B)
	double u_A, u_B, u, u_norm;
	u_A = u_statistic(df, 0.0, n);
	u_B = u_statistic(df, 1.0, n);
	if(u_A < u_B){
		u = u_A;
	} else{
		u = u_B;
	}
	u_norm = u_norm_approx(u, n);
	
	return(0);
}

// fill vector with univariate normal random variables
void generate_normal_vector(double A[], int n, double mu, double var, int seed){

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	int i;
	for(i=0; i < n; i++){
		double x = gsl_ran_gaussian(r, sqrt(var));
		x += mu;
		A[i] = x;
	}

	gsl_rng_free(r);
}

// fill vector with univariate normal random variables
// each variable has a 2% chance of having 10x mu value
void generate_contaminated_normal_vector(double A[], int n, double mu, double var, int seed){

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	int i;
	for(i=0; i < n; i++){

		double x = gsl_ran_gaussian(r, sqrt(var));

		double u = gsl_rng_uniform(r);

		if(u <= 0.02){
			x += (10*mu);
		} else{
			x += mu;
		}
		A[i] = x;
	}

	gsl_rng_free(r);
}

// fill vector with exponential random variables
void generate_exp_vector(double A[], int n, double mu, int seed){

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	int i;
	for(i=0; i < n; i++){
		double x = gsl_ran_exponential(r, mu);
		A[i] = x;
	}

	gsl_rng_free(r);
}

// Welch's t-statistic
// two independent (i.e. non-paired) samples, no assumption of equal variance
// assumes equal sample sizes
double t_statistic(double A[], double B[], int n){

	double mean_A, mean_B, var_A, var_B, t;

	mean_A = gsl_stats_mean(A, 1, n);
	mean_B = gsl_stats_mean(B, 1, n);

	var_A = gsl_stats_variance(A, 1, n);
	var_B = gsl_stats_variance(B, 1, n);

	t = (mean_A - mean_B) / sqrt( (var_A / (double)n) + (var_B / (double)n) );

	return(t);
}

// Welch-Satterthwaite equation to compute degrees of freedom for t-statistic
double t_df(double A[], double B[], int n){

	double var_A, var_B, num, den, df;

	var_A = gsl_stats_variance(A, 1, n);
	var_B = gsl_stats_variance(B, 1, n);
	
	num = pow(((var_A / (double)n) + (var_B / (double)n)), 2);

	den = (pow((var_A / (double)n), 2) / ((double)n - 1)) + (pow((var_B / (double)n), 2) / ((double)n - 1));

	df = num / den;

	return(df);
}

//two-sided p-value for specified t-statistic and df
double t_sig(double t_statistic, double df){

	double t, above, below, two_sided_p;
	
	t = fabs(t_statistic);

	above = gsl_cdf_tdist_Q(t, df);
	below = gsl_cdf_tdist_P(-t, df);
	
	two_sided_p = above + below;

	return(two_sided_p);
}

// copy values from arrays to data structure
void populate_df(struct DF df, double A[], double B[], int n){

	int i;

	for(i = 0; i < n; i++){
		df.Observations[i] = A[i];
		df.Source[i] = 0.0;
	}

	for(i = n; i < (2 * n); i++){
		df.Observations[i] = B[i - n];
		df.Source[i] = 1.0;
	}
}

//bubblesort
//computationally inefficient 
void bsort(struct DF df, int n){

	int i, newn;
	int two_n = 2 * n;
	double templf;
	
	do{
		newn = 0;

		for(i=1; i < two_n; i++){

			if(df.Observations[i-1] > df.Observations[i]){

				//sort Observations in ascending order
				templf = df.Observations[i-1];
				df.Observations[i-1] = df.Observations[i];
				df.Observations[i] = templf;

				//sort Source by Observations
				templf = df.Source[i-1];
				df.Source[i-1] = df.Source[i];
				df.Source[i] = templf;

				newn = i;
			}
		}
		n = newn;
	} while( n != 0);
}

//rank observations in ascending order
void assign_rank(struct DF df, int n){

	int i;
	int two_n = 2 * n;

	for(i=0; i < two_n; i++){
		df.Rank[i] = (double)i + 1;
	}
}

// U-statistic for Wilcoxon rank-sum test
// two independent samples, no parametric assumptions
// assumes equal sample sizes
double u_statistic(struct DF df, double source, int n){

	int two_n = 2 * n;
	int r = 0;
	int i, u;

	for(i=0; i < two_n; i++){
		if(df.Source[i] == source){
			r += df.Rank[i];
		}
	}

	u = r - ((n * (n+1)) / 2);

	return(u);
}

// standardized value of U, approximately normally distributed for sufficiently large n
// assumes equal sample sizes
double u_norm_approx(double u, int n){

	double mu, std, z;

	mu = pow((double)n, 2) / 2;

	std = sqrt((pow((double)n, 2) * ((2 * (double)n) + 1)) / 12);

	z = (u - mu) / std;

	return(z);
}
