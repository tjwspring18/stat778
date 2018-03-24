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

/*****             FUNCTION PROTOTYPES            *******/
struct DF{
	double *Observations; 
	double *Source;
	double *Rank;
	size_t size;
};
struct DF* makeDF(size_t sz);
void deleteDF(struct DF *df);
void generate_normal_vector(double A[], int n, double mu, double var, int seed);
void generate_exp_vector(double A[], int n, double mu, int seed);
void generate_contaminated_normal_vector(double A[], int n, double mu, double var, int seed);
double t_statistic(double A[], double B[], int n);
double t_df(double A[], double B[], int n);
double t_sig(double t_statistic, double df);
int t_test(double A[], double B[], int n, double alpha);
void populate_df(struct DF *df, double A[], double B[], int n);
void bsort(struct DF *df, int n);
void assign_rank(struct DF *df, int n);
double u_statistic(struct DF *df, double source, int n);
double u_norm_approx(double u, int n);
double norm_sig(double zu);
int wrs_test(struct DF *df, int n, double alpha);
void run_simulation(double a, int d, int n, double m1, double s1, double m2, double s2, int seed);

/*****             MAIN            *******/
int main(){

	//set up rng
	int t = time(NULL);
	srand(t);
	int seed;

	//simulation specifications to test:

	//Number in each group
	int N[] = {25, 50, 100};
	
	//Distribution: 1 = normal, 2 = contaminated normal, 3 = exponential
	int D[] = {1,2,3};

	//Group 2 means (group 1 mean always is 1.0)
	//for exponential distribution this also is var
	//for normal distributions, var always is 1.0
	double M2[] = {1.0, 1.05, 1.10, 1.2};

	//print header
	printf("alpha,distribution,n,m1,s1,m2,s2,I_t,II_t,I_w,II_w\n");

	//execute 1000 iterations for every specification
	//prints to STDOUT
	int n, d, m, i;
	for(n=0; n<3; n++){
		for(d=0; d<3; d++){
			for(m=0; m<4; m++){
				for(i=0; i<1000; i++){

					seed = rand();

					if(D[d] != 3){
						run_simulation(0.05, D[d], N[n],
								1.0, 1.0, M2[m],
								1.0, seed);
					} else{
						run_simulation(0.05, D[d], N[n],
								1.0, 1.0, M2[m],
								M2[m], seed);
					}
				}
			}
		}
	}

	return(0);
}


/*****             FUNCTION DEFINITIONS            *******/

// dynamically create instance df of struct type DF
struct DF* makeDF(size_t sz){
	struct DF *df = malloc(sizeof(struct DF));
	df->Observations = malloc(sz * sizeof(double));
	df->Source = malloc(sz * sizeof(double));
	df->Rank = malloc(sz * sizeof(double));
	df->size = sz;
	return(df);
}

// free all memory associated with df
void deleteDF(struct DF *df){
	free(df->Observations);
	free(df->Source);
	free(df->Rank);
	free(df);
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

// t-test
// returns 0 if fail to reject null
// return 1 if reject null
int t_test(double A[], double B[], int n, double alpha){

	double t, df, sig;

	t = t_statistic(A, B, n);
	df = t_df(A, B, n);
	sig = t_sig(t, df);

	if(sig > alpha){
		return(0); //fail to reject null
	} else{ 
		return(1); //reject null
	}
}

// copy values from arrays to data structure
void populate_df(struct DF *df, double A[], double B[], int n){

	int i;

	for(i = 0; i < n; i++){
		df->Observations[i] = A[i];
		df->Source[i] = 0.0;
	}

	for(i = n; i < (2 * n); i++){
		df->Observations[i] = B[i - n];
		df->Source[i] = 1.0;
	}
}

//bubblesort
//computationally inefficient 
void bsort(struct DF *df, int n){

	int i, newn;
	int two_n = 2 * n;
	double templf;
	
	do{
		newn = 0;

		for(i=1; i < two_n; i++){

			if(df->Observations[i-1] > df->Observations[i]){

				//sort Observations in ascending order
				templf = df->Observations[i-1];
				df->Observations[i-1] = df->Observations[i];
				df->Observations[i] = templf;

				//sort Source by Observations
				templf = df->Source[i-1];
				df->Source[i-1] = df->Source[i];
				df->Source[i] = templf;

				newn = i;
			}
		}
		n = newn;
	} while( n != 0);
}

//rank observations in ascending order
void assign_rank(struct DF *df, int n){

	int i;
	int two_n = 2 * n;

	for(i=0; i < two_n; i++){
		df->Rank[i] = (double)i + 1;
	}
}

// U-statistic for Wilcoxon rank-sum test
// two independent samples, no parametric assumptions
// assumes equal sample sizes
double u_statistic(struct DF *df, double source, int n){

	int two_n = 2 * n;
	int r = 0;
	int i, u;

	for(i=0; i < two_n; i++){
		if(df->Source[i] == source){
			r += df->Rank[i];
		}
	}

	u = r - (((double)n * ((double)n+1)) / 2);

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

// two-sided p-value for n(0,1) random variable
double norm_sig(double zu){

	double z, above, below, two_sided_p;
	
	z = fabs(zu);

	above = gsl_cdf_gaussian_Q(z, 1.0);
	below = gsl_cdf_gaussian_P(-z, 1.0);
	
	two_sided_p = above + below;

	return(two_sided_p);
}

// Wilcoxon rank-sum test
// returns 0 if fail to reject null
// returns 1 if reject null
int wrs_test(struct DF *df, int n, double alpha){

	double u, u_norm, sig;

	//sort Observations and Source in ascending order
	bsort(df, n);

	//assign Rank
	assign_rank(df, n);
	
	//calculate u statistic
	u = u_statistic(df, 0.0, n);

	//take normal approximation
	u_norm = u_norm_approx(u, n);

	//two-sided p-value
	sig = norm_sig(u_norm);

	if(sig > alpha){
		return(0); //fail to reject null
	} else{
		return(1); //reject null
	}
}

// a = alpha
// d = distribution. 1 = normal, 2 = contaminated normal, 3 = exponential
// n = number of samples in each group (implicitly assumes equal group size)
// m1 = mean of group 1
// s1 = variance of group 1. not used if d = 3.
// m2 = mean of group 2.
// s2 = variance of group 2. not used if d = 3.
void run_simulation(double a, int d, int n, double m1, double s1, double m2, double s2, int seed){

	int t, w, I_t, II_t, I_w, II_w;

	//generate vectors in which to store RVs
	//we will use these for t-test
	double *A = malloc(n * sizeof(double));
	double *B = malloc(n * sizeof(double));


	// generate random variables, store in A and B 
	if(d == 1){
		generate_normal_vector(A, n, m1, s1, seed);
		generate_normal_vector(B, n, m2, s2, seed);
	} else{
		if(d == 2){
			generate_contaminated_normal_vector(A, n, m1, s1, seed);
			generate_contaminated_normal_vector(B, n, m1, s1, seed);
		} else{
			generate_exp_vector(A, n, m1, seed);
			generate_exp_vector(B, n, m1, seed);
		}
	}

	// initialize data structure
	// copy data from A and B to data structure
	// we will use this for Wilcoxon rank-sum test
	struct DF *df;
	df = makeDF(2*n);
	populate_df(df, A, B, n);

	//conduct t-test
	t = t_test(A, B, n, a);

	//conduct Wilcoxon rank-sum test
	w = wrs_test(df, n, a);

	//check to see if type I error committed
	//by t-test
	if((m1 == m2) & (t == 1)){
		I_t = 1;
	} else{
		I_t = 0;
	}
	//and by Wilcoxon rank sum test
	if((m1 == m2) & (w == 1)){
		I_w = 1;
	} else{
		I_w = 0;
	}

	//check to see if type II error committed 
	//by t-test
	if( (m1 != m2) & (t == 0)){
		II_t = 1;
	} else{
		II_t = 0;
	}

	//and by Wilcoxon rank sum test
	if( (m1 != m2) & (w == 0)){
		II_w = 1;
	} else{
		II_w = 0;
	}

	//print simulation parameters and results
	printf("%lf,%d,%d,%lf,%lf,%lf,%lf,%d,%d,%d,%d\n",
			a, d, n, m1, s1, m2, s2, I_t, II_t, I_w, II_w);

	// free memory
	free(A);
	free(B);
	deleteDF(df);
}
