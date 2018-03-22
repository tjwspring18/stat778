#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_statistics_double.h>
#include<gsl/gsl_cdf.h>

// Function prototypes and struct definitions
void generate_normal_vector(double A[], int n, double mu, double var, int seed);
void generate_exp_vector(double A[], int n, double mu, int seed);
void generate_contaminated_normal_vector(double A[], int n, double mu, double var, int seed);
double t_statistic(double A[], double B[], int n_A, int n_B);
double t_df(double A[], double B[], int n_A, int n_B);
double t_sig(double t_statistic, double df);
struct DF{
	double *Observations; 
	double *Source;
	double *Rank;
};
/* All of the below are TODO */
void populate_df(struct DF df, double A[], double B[], int n_A, int n_B);
void bsort(struct DF df, int n);
void assign_rank(struct DF df, int n);
double u_statistic(struct DF df, double source, int n);

/*
   Wilcoxon rank sum test

   Make 2n x 3 array

   Put RVs from A and B in first column

   Put marker of whether came from A and B in second column

   Sort array

   Put rank of observation in 3rd column

   Add up sum of ranks for A: call this U1

   Use smaller of U1 or U2

   Use normal approximation?

 */

int main(){

	//number in each group
	int n = 10;

	//initialize random number generator
	time_t t;
	srand((unsigned) time(&t));


	double *A = malloc(n * sizeof(double));
	double *B = malloc(n * sizeof(double));
	generate_normal_vector(A, n, 0, 1, rand());
	generate_normal_vector(B, n, 0, 1, rand());
	
	/*
	double t_score; 
	t_score = t_statistic(A, B, n, n);
	printf("%lf\n", t_score);
	*/

	//something wonky here
	/*
	double val;
	val= t_test(t_score, n);
	printf("%lf\n", val);
	*/

	/* Begin testing */
	/*
	for(i=0; i < 20; i++){
		printf("%lf\n", A[i]);
	}

	double mean = gsl_stats_mean(A, 1, n);
	printf("%lf\n", mean);
	*/

	return(0);
}

void generate_normal_vector(double A[], int n, double mu, double var, int seed){

	// set up rng
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	// generate normal RVs and stash in A[]
	int i;
	for(i=0; i < n; i++){
		double x = gsl_ran_gaussian(r, sqrt(var));
		x += mu;
		A[i] = x;
	}

	gsl_rng_free(r);
}

void generate_contaminated_normal_vector(double A[], int n, double mu, double var, int seed){

	// set up rng
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	// generate normal RVs and stash in A[]
	// 2% chance of abnormally large value
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

void generate_exp_vector(double A[], int n, double mu, int seed){

	// set up rng
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	// generate exponential RVs and stash in A[]
	int i;
	for(i=0; i < n; i++){
		double x = gsl_ran_exponential(r, mu);
		A[i] = x;
	}

	gsl_rng_free(r);
}

// Calculate t-statistic by Welch's method
double t_statistic(double A[], double B[], int n_A, int n_B){

	double mean_A = gsl_stats_mean(A, 1, n_A);
	double mean_B = gsl_stats_mean(B, 1, n_B);

	double var_A = gsl_stats_variance_m(A, 1, n_A, mean_A);
	double var_B = gsl_stats_variance_m(B, 1, n_B, mean_A);

	double t = (mean_A - mean_B) / sqrt( (var_A / (double)n_A) + (var_B / (double)n_B) );

	return(t);
}

// Calculate degrees of freedom by Welch-Satterthwaite equation
double t_df(double A[], double B[], int n_A, int n_B){

	double var_A = gsl_stats_variance(A, 1, n_A);

	double var_B = gsl_stats_variance(B, 1, n_B);
	
	double num = pow(((var_A / (double)n_A) + (var_B / (double)n_B)), 2)

	double den = (pow((var_A / (double)n_A), 2) / ((double)n_A - 1)) + (pow((var_B / (double)n_B), 2) / ((double)n_B - 1))

	double df = num / den;

	return(df)
}

// Think about this more
// Needs to be two-sided
double t_sig(double t_statistic, double df){
	double t = fabs(t_statistic);
	double x;
	x = gsl_cdf_tdist_Q(t, df);
	return(x);
}
