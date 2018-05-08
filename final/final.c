#include<stdio.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_statistics_double.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>

/*
   Tom Wallace <twalla11@masonlive.gmu.edu>
   STAT 778
   Spring 2018
   Final Project
 */

/************************************************************************
                  DATA STRUCTURES AND FUNCTION PROTOTYPES
 ***********************************************************************/
struct Data{
	double *t;     // follow-up time
	double *y;     // count variable
	double *x;     // control or treatment group
	size_t size;
};
struct Data* makeData(size_t sz);
void deleteData(struct Data *data);
double min(double a, double b);
void assignT(struct Data *data, double lambda_d);
void assignX(struct Data *data, double group);
void assignY(struct Data *data, double p, double rr);
void runSimulation(int n_assigned, double lambda_d, double p_1, double rr_1, double p_2, double rr_2);
double nr(struct Data *control, struct Data *treatment);
double f_1(struct Data *control, struct Data *treatment, double b_0, double b_1);
double f_2(struct Data *control, struct Data *treatment, double b_0, double b_1);
double f_1_b0(struct Data *control, struct Data *treatment, double b_0, double b_1);
double f_1_b1(struct Data *control, struct Data *treatment, double b_0, double b_1);
double f_2_b1(struct Data *control, struct Data *treatment, double b_0, double b_1);

/************************************************************************
                                 MAIN
 ************************************************************************/
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
	data->t = malloc(sz * sizeof(double));
	data->y = malloc(sz * sizeof(double));
	data->x = malloc(sz * sizeof(double));
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
		data->t[i] = min((double)k, 1.0);
	}

	gsl_rng_free(r);
}

double min(double a, double b){
	if(a < b){
		return(a);
	} else{
		return(b);
	}
}

void assignX(struct Data *data, double group){
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

	double b_1;

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
	/*
	double mean = gsl_stats_mean(control->y, 1, control->size);
	double var = gsl_stats_variance(control->y, 1, control->size);
	printf("%lf %lf\n", mean, var);

	mean = gsl_stats_mean(treatment->y, 1, treatment->size);
	var = gsl_stats_variance(treatment->y, 1, treatment->size);
	printf("%lf %lf\n", mean, var);
	*/

	//functions to calculate hat_beta
	//TODO: this ain't working
	b_1 = nr(control, treatment);
	printf("%lf\n", b_1);
	
	//functions to calculate hat_w0

	//do test

	deleteData(control);
	deleteData(treatment);
}

double nr(struct Data *control, struct Data *treatment){

	int s;
	double e = 1;
	double b_0, b_1;

	//allocate matrices we'll need for Newton-Raphson
	gsl_matrix *B = gsl_matrix_alloc(2,1);
	gsl_matrix *Bnew = gsl_matrix_alloc(2,1);
	gsl_matrix *J = gsl_matrix_alloc(2,1);
	gsl_matrix *H = gsl_matrix_alloc(2,2);
	gsl_matrix *H_inv = gsl_matrix_alloc(2,2);
	gsl_matrix *prod = gsl_matrix_alloc(2,1);
	gsl_matrix *diff = gsl_matrix_alloc(2,1);

	//initial guesses for betas are 0
	gsl_matrix_set(B, 0, 0, 0);
	gsl_matrix_set(B, 1, 0, 0);

	//iterate
	while(e > 10e-6){

		//get current values of beta
		b_0 = gsl_matrix_get(B, 0, 0);
		b_1 = gsl_matrix_get(B, 1, 0);
		

		//evaluate functions at current estimates
		gsl_matrix_set(J, 0, 0, f_1(control, treatment, b_0, b_1));
		gsl_matrix_set(J, 1, 0, f_2(control, treatment, b_0, b_1));
		printf("%lf\n", gsl_matrix_get(J, 0, 0));
		printf("%lf\n", gsl_matrix_get(J, 1, 0));

		//evaluate derivatives at current estimates 
		gsl_matrix_set(H, 0, 0, f_1_b0(control, treatment, b_0, b_1));
		gsl_matrix_set(H, 0, 1, f_1_b1(control, treatment, b_0, b_1));
		gsl_matrix_set(H, 1, 0, f_1_b1(control, treatment, b_0, b_1));
		gsl_matrix_set(H, 1, 1, f_2_b1(control, treatment, b_0, b_1));

		/*
		printf("%lf\n", gsl_matrix_get(H, 0, 0));
		printf("%lf\n", gsl_matrix_get(H, 0, 1));
		printf("%lf\n", gsl_matrix_get(H, 1, 0));
		printf("%lf\n", gsl_matrix_get(H, 1, 1));
		*/

		//take inverse of H via LU decomposition
		gsl_permutation *p = gsl_permutation_alloc(2);
		gsl_linalg_LU_decomp(H, p, &s);
		gsl_linalg_LU_invert(H, p, H_inv);
		gsl_permutation_free(p);

		//multiply H_inv and J
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, H_inv, J, 0.0, prod);

		//new betas = current betas - prod
		gsl_matrix_set(Bnew, 0, 0, (gsl_matrix_get(B, 0, 0) - gsl_matrix_get(prod, 0, 0)));
		gsl_matrix_set(Bnew, 1, 0, (gsl_matrix_get(B, 1, 0) - gsl_matrix_get(prod, 1, 0)));

		//calculate difference between B and Bnew
		gsl_matrix_set(diff, 0, 0, (gsl_matrix_get(B, 0, 0) - gsl_matrix_get(Bnew, 0, 0)));
		gsl_matrix_set(diff, 1, 0, (gsl_matrix_get(B, 1, 0) - gsl_matrix_get(Bnew, 1, 0)));

		//set beta equal to new beta
		gsl_matrix_sub(B, prod);

		//get error
		e = fabs(gsl_matrix_max(diff));
	}

	//free memory
	gsl_matrix_free(B);
	gsl_matrix_free(Bnew);
	gsl_matrix_free(J);
	gsl_matrix_free(H);
	gsl_matrix_free(H_inv);
	gsl_matrix_free(prod);
	gsl_matrix_free(diff);

	return(b_1);
}

double f_1(struct Data *control, struct Data *treatment, double b_0, double b_1){
	int i;
	int n_1 = control->size;
	int n_2 = treatment->size;

	double q_1, q_2 = 0;
	double q;

	for(i = 0; i < n_1; i++){
		q_1 += (control->y[i] - (control->t[i] * exp(b_0 + (b_1 * control->x[i]))));
	}

	for(i = 0; i < n_2; i++){
		q_2 += (treatment->y[i] - (treatment->t[i] * exp(b_0 + (b_1 * treatment->x[i]))));
	}

	q = q_1 + q_2;

	return(q);
}

double f_2(struct Data *control, struct Data *treatment, double b_0, double b_1){

	int i;
	int n_1 = control->size;
	int n_2 = treatment->size;

	double q_1, q_2 = 0;
	double q;

	for(i = 0; i < n_1; i++){
		q_1 += (control->x[i]*(control->y[i] - (control->t[i] * exp(b_0 + (b_1 * control->x[i])))));
	}

	for(i = 0; i < n_2; i++){
		q_2 += (treatment->x[i] * (treatment->y[i] - (treatment->t[i] * exp(b_0 + (b_1 * treatment->x[i])))));
	}

	q = q_1 + q_2;

	return(q);
}

double f_1_b0(struct Data *control, struct Data *treatment, double b_0, double b_1){

	int i;
	int n_1 = control->size;
	int n_2 = treatment->size;

	double q_1, q_2 = 0;
	double q;

	for(i = 0; i < n_1; i++){
		q_1 += (control->t[i] * exp(b_0 + (b_1 * control->x[i])));
	}

	for(i = 0; i < n_2; i++){
		q_2 += (treatment->t[i] * exp(b_0 + (b_1 * treatment->x[i])));
	}

	q = -1.0 * (q_1 + q_2);

	return(q);
}

double f_1_b1(struct Data *control, struct Data *treatment, double b_0, double b_1){

	int i;
	int n_1 = control->size;
	int n_2 = treatment->size;

	double q_1, q_2 = 0;
	double q;

	for(i = 0; i < n_1; i++){
		q_1 += (control->t[i] * control->x[i] * exp(b_0 + (b_1 * control->x[i])));
	}

	for(i = 0; i < n_2; i++){
		q_2 += (treatment->t[i] * treatment->x[i] * exp(b_0 + (b_1 * treatment->x[i])));
	}

	q = -1.0 * (q_1 + q_2);

	return(q);
}

double f_2_b1(struct Data *control, struct Data *treatment, double b_0, double b_1){

	int i;
	int n_1 = control->size;
	int n_2 = treatment->size;

	double q_1, q_2 = 0;
	double q;

	for(i = 0; i < n_1; i++){
		q_1 += (control->t[i] * pow(control->x[i], 2) * exp(b_0 + (b_1 * control->x[i])));
	}

	for(i = 0; i < n_2; i++){
		q_2 += (treatment->t[i] * pow(treatment->x[i], 2) * exp(b_0 + (b_1 * treatment->x[i])));
	}

	q = -1.0 * (q_1 + q_2);

	return(q);
}

double phi_tilde_hat(struct Data *control, struct Data *treatment){
	
	int i;
	int n = control->size;
	
	for(i = 0; i < n; i++){
		
	}

}
