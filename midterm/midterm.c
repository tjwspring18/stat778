#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_statistics_double.h>
#include<gsl/gsl_cdf.h>

void generate_normal_vector(double A[], int n, double mu, double var, int seed);
void generate_exp_vector(double A[], int n, double mu, int seed);
void generate_contaminated_normal_vector(double A[], int n, double mu, double
		var, int seed);
double t_statistic(double A[], double B[], int n_A, int n_B);
double t_test(double t_statistic, int n);
/*
TODO: fix t_test
TODO: Wilxocon
*/
int main(){


	//initialize random number generator
	time_t t;
	srand((unsigned) time(&t));

	int n = 10;

	double *A = malloc(n * sizeof(double));
	double *B = malloc(n * sizeof(double));
	generate_normal_vector(A, n, 0, 1, rand());
	generate_normal_vector(B, n, 0, 1, rand());
	
	double t_score; 
	t_score = t_statistic(A, B, n, n);
	printf("%lf\n", t_score);

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

	// generate RVs and stash in A[]
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

		if(u < 0.02){
			x += 10*mu;
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

double t_statistic(double A[], double B[], int n_A, int n_B){

	double mean_A = gsl_stats_mean(A, 1, n_A);
	double mean_B = gsl_stats_mean(B, 1, n_B);

	double var_A = gsl_stats_variance_m(A, 1, n_A, mean_A);
	double var_B = gsl_stats_variance_m(B, 1, n_B, mean_A);

	double t = (mean_A - mean_B) / sqrt( (var_A / n_A) + (var_B / n_B) );

	return(t);
}

double t_test(double t_statistic, int n){
	double t = abs(t_statistic);
	double x;
	x = gsl_cdf_tdist_Q(t, n-1);
	return(x);
}
