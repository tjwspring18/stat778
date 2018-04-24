#include<stdio.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_statistics_double.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

/*
   Tom Wallace
   STAT 778
   Final Project
 */

/************************************************************************
                  DATA STRUCTURES AND FUNCTION PROTOTYPES
 ***********************************************************************/
struct Data{
	int *t;     // follow-up time
	double *y;     // count variable
	int *x;     // control or treatment group
	size_t size;
};
struct Data* makeData(size_t sz);
void deleteData(struct Data *data);
int min(int a, int b);
void assignT(struct Data *data, double lambda_d);
void assignX(struct Data *data, int group);
void assignY(struct Data *data, double p, double rr);
void runSimulation(int n_assigned, double lambda_d, double p_1, double rr_1, double p_2, double
		rr_2);

/************************************************************************
                                 MAIN
 ***********************************************************************/
int main(void){

	int n_assigned = 1000;
	double lambda_d = 1 / 0.356;
	double p_1 = 1-0.46352;
	double rr_1 = 1.44625;
	double p_2 = 1 - 0.49239;
	double rr_2 = 1.0309;

	runSimulation(n_assigned, lambda_d, p_1, rr_1, p_2, rr_2);

}

/************************************************************************
                            FUNCTION DEFINITIONS
 ***********************************************************************/

struct Data* makeData(size_t sz){
	struct Data *data = malloc(sizeof(struct Data));
	data->t = malloc(sz * sizeof(int));
	data->y = malloc(sz * sizeof(double));
	data->x = malloc(sz * sizeof(int));
	data->size = sz;
	return(data);
}

void deleteData(struct Data *data){
	free(data->t);
	free(data->x);
	free(data->y);
	free(data);
}

void assignT(struct Data *data, double lambda_d){

	int i;
	int n = data->size;

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());

	for(i = 0; i < n; i++){
		int k = gsl_ran_exponential(r, lambda_d);
		data->t[i] = min(k, 1);
	}

	gsl_rng_free(r);
}

int min(int a, int b){
	if(a < b){
		return(a);
	} else{
		return(b);
	}
}

void assignX(struct Data *data, int group){
	int n = data->size;
	int i;
	for(i = 0; i < n; i++){
		data->x[i] = group;
	}
}

void assignY(struct Data *data, double p, double rr){

	int n = data->size;
	int i, c;

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());

	for(i = 0; i < n; i++){
		c = gsl_ran_negative_binomial(r, p, rr);
		data->y[i] = (double)c;
	}

	gsl_rng_free(r);
}

void runSimulation(int n_assigned, double lambda_d, double p_1, double rr_1, double p_2, double rr_2){

	srand(time(NULL));

	struct Data *control;
	struct Data *treatment;

	// generate two groups
	control = makeData(n_assigned);
	treatment = makeData(n_assigned);

	// covariate X is indicator of group membership
	assignX(control, 0);
	assignX(treatment, 1);

	// follow up time follows exponential distribution
	// about 30% of subjects drop out early
	assignT(control, lambda_d);
	assignT(treatment, lambda_d);

	// outcome variable (count) Y
	assignY(control, p_1, rr_1);
	assignY(treatment, p_2, rr_2);

	/* Replace p with 1-p */
	double mean = gsl_stats_mean(control->y, 1, control->size);
	double var = gsl_stats_variance(control->y, 1, control->size);
	printf("%lf %lf\n", mean, var);

	mean = gsl_stats_mean(treatment->y, 1, treatment->size);
	var = gsl_stats_variance(treatment->y, 1, treatment->size);
	printf("%lf %lf\n", mean, var);

	//functions to calculate hat_beta

	//functions to calculate hat_w0

	//do test

	deleteData(control);
	deleteData(treatment);
}
