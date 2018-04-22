#include<stdio.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

int min(int a, int b){
	if(a < b){
		return(a);
	} else{
		return(b);
	}
}
/*
   row 1: p = 0.46352, r = 1.44
 */

int main(void){
	srand(time(NULL));

	int i;
	double lambda_d = 1 / 0.356;

	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, rand());

	for(i=0; i<100; i++){
		int k = gsl_ran_exponential(r, lambda_d);
		printf("%d\n", k);
	}

	gsl_rng_free(r);

}

