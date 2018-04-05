#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#define ABORT_MESSAGE "Program aborting"
#define BAD_ARGS "Wrong number of arguments passed"
#define FNF "File not found"

/*
   Tom Wallace <twalla11@masonlive.gmu.edu>
   STAT 778
   Spring 2018
   HW 3
 */

/**********************************************************************
  DATA STRUCTURES AND FUNCTION PROTOTYPES
 **********************************************************************/
struct Data{
	double *t;  // event time
	int *e;     // event indicator
	double *x1; // covariate 1
	double *x2; // covariate 2
	size_t size;
};
struct Data* makeData(size_t sz);
void deleteData(struct Data *data);
int count_lines(char *f); 
void read_data(char *f, struct Data *data);
void bsort(struct Data *data, int n);
double deriv_first_x1(struct Data *data, double b1, double b2);
double deriv_first_x2(struct Data *data, double b1, double b2);
double deriv_second_x1_x2(struct Data *data, double b1, double b2);
double deriv_second_x1_x1(struct Data *data, double b1, double b2);
double deriv_second_x2_x2(struct Data *data, double b1, double b2);

/**********************************************************************
  MAIN
 **********************************************************************/
int main(int argc, char *argv[]){

	static char *in;
	static int n;
	int i, s;
	struct Data *data;
	double e = 1;
	double b1, b2;

	// make sure input file argument passed
	if(argc != 2){
		printf("%s\n", BAD_ARGS);
		printf("%s\n", ABORT_MESSAGE);
		exit(1);
	} else{
		//get name of input file
		in = argv[1];

		//check that it exists
		if(access(in, F_OK) == -1){
			printf("%s\n", FNF);
			printf("%s\n", ABORT_MESSAGE);
			exit(1);
		} else {
			printf("Reading file %s\n", in);

			//count number of observations
			n = count_lines(in);

			//dynamically create data structure
			data = makeData(n);

			//read data from file into data structure
			read_data(in, data);
			printf("Read %d observations\n", n);

			//sort data by time
			bsort(data, n);

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

			//Newton Raphson 
			while(e > 10e-6){

				//get current values of beta
				b1 = gsl_matrix_get(B, 0, 0);
				b2 = gsl_matrix_get(B, 1, 0);

				//evaluate first derivatives at current
				//estimates of beta1 and beta2
				gsl_matrix_set(J, 0, 0, deriv_first_x1(data, b1, b2));
				gsl_matrix_set(J, 1, 0, deriv_first_x2(data, b1, b2));

				//evaluate second derivatives at current
				//estimates of beta1 and beta2
				gsl_matrix_set(H, 0, 0, deriv_second_x1_x1(data, b1, b2));
				gsl_matrix_set(H, 0, 1, deriv_second_x1_x2(data, b1, b2));
				gsl_matrix_set(H, 1, 0, deriv_second_x1_x2(data, b1, b2));
				gsl_matrix_set(H, 1, 1, deriv_second_x2_x2(data, b1, b2));

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

			/*
			   then get SE
			   evaluate Hessian at final estimates of beta-hat
			   ''the square roots of the diagonal elements of the
			   inverse of the Hessian are are the estimated standard
			   errors''

			 */
			//free matrices
			gsl_matrix_free(B);
			gsl_matrix_free(Bnew);
			gsl_matrix_free(J);
			gsl_matrix_free(H);
			gsl_matrix_free(H_inv);
			gsl_matrix_free(prod);

			//delete Data
			deleteData(data);
			printf("%lf %lf\n", b1, b2);

			return(0);
		}
	}
}

/**********************************************************************
  FUNCTION DEFINITIONS
 **********************************************************************/

// dynamically create instance data of struct type Data
struct Data* makeData(size_t sz){
	struct Data *data = malloc(sizeof(struct Data));
	data->t = malloc(sz * sizeof(double));
	data->e = malloc(sz * sizeof(int));
	data->x1 = malloc(sz * sizeof(double));
	data->x2 = malloc(sz * sizeof(double));
	data->size = sz;
	return(data);
}

// delete instance of Data
void deleteData(struct Data *data){
	free(data->t);
	free(data->e);
	free(data->x1);
	free(data->x2);
	free(data);
}

// count number of lines in input file
int count_lines(char *f){

	int lines = 0;
	char ch;
	FILE *fp = fopen(f, "r");

	while(!feof(fp)){
		ch = fgetc(fp);
		if(ch == '\n'){
			lines++;
		}
	}

	fclose(fp);
	return(lines);
}

//read data from input file into data structure
void read_data(char *f, struct Data *data){

	FILE *fp = fopen(f, "r");

	double t;
	int e;
	double x1;
	double x2;

	int i = 0;
	while(fscanf(fp, "%lf %d %lf %lf", &t, &e, &x1, &x2) != EOF){
		data->t[i] = t;
		data->e[i] = e;
		data->x1[i] = x1;
		data->x2[i] = x2;
		i++;
	}
	fclose(fp);
}

//bubblesort
//computationally inefficient 
void bsort(struct Data *data, int n){

	int i, newn;
	double templf;

	do{
		newn = 0;

		for(i=1; i < n; i++){

			if(data->t[i-1] > data->t[i]){

				//sort t in ascending order
				templf = data->t[i-1];
				data->t[i-1] = data->t[i];
				data->t[i] = templf;

				//sort e by t
				templf = data->e[i-1];
				data->e[i-1] = data->e[i];
				data->e[i] = templf;

				//sort x1 by t
				templf = data->x1[i-1];
				data->x1[i-1] = data->x1[i];
				data->x1[i] = templf;

				//sort x2 by t
				templf = data->x2[i-1];
				data->x2[i-1] = data->x2[i];
				data->x2[i] = templf;

				newn = i;

			}
		}
		n = newn;
	} while( n != 0);
}

double deriv_first_x1(struct Data *data, double b1, double b2){

	int i, j, n;
	double d = 0;
	double num, den;

	//get number of observations
	n = (int)data->size;

	for(i=0; i < n; i++){

		num = 0;
		den = 0;

		//if observation is not censored
		if(data->e[i] == 1){

			d += data->x1[i];

			for(j=i; j < n; j++){
				num += (data->x1[j] * exp((b1 * data->x1[j]) + (b2 * data->x2[j])));
				den += exp((b1 * data->x1[j]) + (b2 * data->x2[j]));
			}

			d -= (num / den);

		} 
	}

	return(d);
}

double deriv_first_x2(struct Data *data, double b1, double b2){

	int i, j, n;
	double d = 0;
	double num, den;

	//get number of observations
	n = (int)data->size;

	for(i=0; i < n; i++){

		num = 0;
		den = 0;

		//if observation is not censored
		if(data->e[i] == 1){

			d += data->x2[i];

			for(j=i; j < n; j++){
				num += (data->x2[j] * exp((b1 * data->x1[j]) + (b2 * data->x2[j])));
				den += exp((b1 * data->x1[j]) + (b2 * data->x2[j]));
			}

			d -= (num / den);

		} 
	}

	return(d);
}

double deriv_second_x1_x2(struct Data *data, double b1, double b2){

	int i, j, n;
	double d = 0;
	double num1, den1, num2a, num2b, den2, e;

	//get number of observations
	n = (int)data->size;

	for(i=0; i<n; i++){

		num1 = 0;
		den1 = 0;
		num2a = 0; 
		num2b = 0;
		den2 = 0;

		//if observation is not censored
		if(data->e[i] == 1){

			for(j=i; j<n; j++){

				e = exp((b1 * data->x1[j]) + (b2 * data->x2[j]));

				num1 += (data->x1[j] * data->x2[j] * e);

				den1 += e;

				num2a += (data->x1[j] * e);

				num2b += (data->x2[j] * e);

				den2 += pow(e, 2);

			}


			d += ((num1 / den1) - ((num2a * num2b) / den2));
		}
	}

	d = (d * -1.0);

	return(d);

}

double deriv_second_x1_x1(struct Data *data, double b1, double b2){

	int i, j, n;
	double d = 0;
	double num1, den1, num2a, num2b, den2, e;

	//get number of observations
	n = (int)data->size;

	for(i=0; i<n; i++){

		num1 = 0;
		den1 = 0;
		num2a = 0; 
		num2b = 0;
		den2 = 0;

		//if observation is not censored
		if(data->e[i] == 1){

			for(j=i; j<n; j++){

				e = exp((b1 * data->x1[j]) + (b2 * data->x2[j]));

				num1 += (data->x1[j] * data->x1[j] * e);

				den1 += e;

				num2a += (data->x1[j] * e);

				num2b += (data->x2[j] * e);

				den2 += pow(e, 2);

			}


			d += ((num1 / den1) - ((num2a * num2b) / den2));
		}
	}

	d = (d * -1.0);

	return(d);

}

double deriv_second_x2_x2(struct Data *data, double b1, double b2){

	int i, j, n;
	double d = 0;
	double num1, den1, num2a, num2b, den2, e;

	//get number of observations
	n = (int)data->size;

	for(i=0; i<n; i++){

		num1 = 0;
		den1 = 0;
		num2a = 0; 
		num2b = 0;
		den2 = 0;

		//if observation is not censored
		if(data->e[i] == 1){

			for(j=i; j<n; j++){

				e = exp((b1 * data->x1[j]) + (b2 * data->x2[j]));

				num1 += (data->x2[j] * data->x2[j] * e);

				den1 += e;

				num2a += (data->x1[j] * e);

				num2b += (data->x2[j] * e);

				den2 += pow(e, 2);

			}


			d += ((num1 / den1) - ((num2a * num2b) / den2));
		}
	}

	d = (d * -1.0);

	return(d);

}
