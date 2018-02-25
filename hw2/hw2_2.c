#include<stdio.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_statistics_double.h>

/*
   gcc -Wall -I/usr/include -c hw2_2.c && gcc -L/usr/lib/i386-linux-gnu hw2_2.o -o hw2_2
*/

int main(){
	//process command line args

	// set up rng
	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	// parameters for normal random variables
	double mu = -0.5;
	double ss = 2.0;

	// number of RVs to generate
	int n=50;

	// create array in which to store normal RVs
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

	// SE of mean
	double se_mean = sqrt(var / n);

	// SE of ss
	double se_var = sqrt((2*pow(var, 2)) / (n-1));

	// 95% CI of mean
	double mean_u = mean + 1.96*se_mean;
	double mean_l = mean - 1.96*se_mean;

	// 95% CI of sss
	double var_u = var + 1.96*se_var;
	double var_l = var - 1.96*se_var;

	/*
	printf("Sample mean: %lf\n", mean);
	printf("Sample variance: %lf\n", var);
	printf("SE(sample mean): %lf\n", se_mean);
	printf("SE(sample variance): %lf\n", se_var);
	*/
}
