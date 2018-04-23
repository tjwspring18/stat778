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
int dropouts(int n_assigned, double lambda_d);
int min(int a, int b);
void assignGroup(struct Data *data, int group);
void assignCount(struct Data *data, double p, double rr);

/************************************************************************
                                 MAIN
 ***********************************************************************/
int main(void){

	srand(time(NULL));

	struct Data *control;
	struct Data *treatment;

	int i;
	int n_observed_control;
	int n_observed_treatment;
	double lambda_d = 1 / 0.356;

	int n_assigned = 1000;

	n_observed_control = dropouts(n_assigned, lambda_d);
	n_observed_treatment = dropouts(n_assigned, lambda_d);

	control = makeData(n_observed_control);
	treatment = makeData(n_observed_treatment);

	assignGroup(control, 0);
	assignGroup(treatment, 0);

	assignCount(control, 1-0.46352, 1.44625);
	assignCount(treatment, 0.5, 1.44625);

/*
   row 1: p = 0.46352, r = 1.44
 */
	/* Replace p with 1-p */
	double mean = gsl_stats_mean(control->y, 1, control->size);
	double var = gsl_stats_variance(control->y, 1, control->size);
	printf("%lf %lf\n", mean, var);

	//functions to calculate hat_beta

	//functions to calculate hat_w0

	//do test

	//record whether type II error occurs
	/*
	for(i = 0; i < data->size; i++){
		printf("%d %d\n", data->x[i], data->y[i]);
	}
	*/
	deleteData(control);
	deleteData(treatment);

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

int dropouts(int n_assigned, double lambda_d){

	int n_observed = 0;
	int i;

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());

	for(i = 0; i < n_assigned; i++){
		int k = gsl_ran_exponential(r, lambda_d);
		k = min(k, 1);
		if(k > 0){
			n_observed++;
		}
	}

	gsl_rng_free(r);

	return(n_observed);
}

int min(int a, int b){
	if(a < b){
		return(a);
	} else{
		return(b);
	}
}

void assignGroup(struct Data *data, int group){
	int n = data->size;
	int i;
	for(i = 0; i < n; i++){
		data->x[i] = group;
	}
}

void assignCount(struct Data *data, double p, double rr){

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
