#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
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

/**********************************************************************
                                MAIN
 **********************************************************************/
int main(int argc, char *argv[]){

	static char *in;
	static int n;
	int i;
	struct Data *data;

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

			double d1_x1, d1_x2, d2_x1_x2;

			d1_x1 = deriv_first_x1(data, 0, 0);

			d1_x2 = deriv_first_x2(data, 0, 0);

			d2_x1_x2 = deriv_second_x1_x2(data, 0, 0);

			printf("%lf %lf %lf\n", d1_x1, d1_x2, d2_x1_x2);

			//deriv_first(data, 1);
			//Newton-Raphson optimization to estimate betas

			/*
			   double J[2];
			   double H[2][2];
			   double b[2]={0,0};
			   double bnew[2]={0,0};
			   double q[2]={0,0}
			   double e1, e2;

			   do{

			   J[0] = deriv_first("x1", data, beta);
			   J[1] = deriv_first("x2", data, beta);

			   H[0][0] = deriv_second(x1, x1, data, beta);
			   H[0][1] = deriv_second(x1, x2, data, beta);
			   H[1][0] = H[0][1];
			   H[1][1] = deriv_second(x2, x2, data, beta);

			   check for all positive eigenvalues

			   take Cholesky composition, get inverse

			   matrix multiply H_inv and J = q

			   bnew[0] = b[0] - q[0]
			   bnew[1] = b[1] - q[1]

			   e1 = bnew[0] - b[0]
			   e2 = bnew[1] - b[1]

			   b[0] = bnew[0]
			   b[1] = bnew[1]
			   
			   } while(e1 > 0.0001 & e2 > 0.0001)

			   then get SE
			   evaluate Hessian at final estimates of beta-hat
			   ''the square roots of the diagonal elements of the
			   inverse of the Hessian are are the estimated standard
			   errors''

			 */

			/*
			for(i=0; i<n; i++){
				printf("%lf %d %lf %lf %lf\n",
						data->t[i],
						data->e[i],
						data->x1[i],
						data->x2[i]);
			}
			*/

			//delete Data
			deleteData(data);

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

				printf("%lf\n", e);

				num1 += (data->x1[j] * data->x2[j] * e);

				den1 += e;

				num2a += (data->x1[j] * e);

				num2b += (data->x2[j] * e);

				num2b += pow(e, 2);

			}


			d += ((num1 / den1) - ((num2a * num2b) / den2));
		}
	}

	d = (d * -1.0);

	return(d);

}
