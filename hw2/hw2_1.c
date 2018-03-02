#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
/*
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>
 */
#define ABORT_MESSAGE "Program aborting"
#define BAD_ARGS "Wrong number of arguments passed"
#define FNF "File not found"

/*
   struct rparams{
   double a;
   double b;
   };

   int rosenbrock_f(const gsl_vector *x, void *params, gsl_vector *f){
   double a = ((struct rparams *) params)->a;
   double b = ((struct rparams *) params)->b;

   const double x0 = gsl_vector_get(x,0);
   const double x1 = gsl_vector_get(x,1);

   const double y0 = a * (1-x0);
   const double y1 = b * (x1 - x0 * x0);

   gsl_vector_set(f, 0, y0);
   gsl_vector_set(f, 1, y1);

   return GSL_SUCCESS;
   }

   int print_state(size_t iter, gsl_multiroot_fsolver *s){
   printf("iter = %3u x = % .3f % .3f "
   "f(x) = %.3e %.3e\n",
   iter,
   gsl_vector_get(s->x, 0),
   gsl_vector_get(s->x, 1),
   gsl_vector_get(s->f, 0),
   gsl_vector_get(s->f, 1));
   }
   int main(void){
   const gsl_multiroot_fsolver_type *T;
   gsl_multiroot_fsolver *s;

   int status;
   size_t i, iter=0;

   const size_t n = 2;
   struct rparams p = {1.0, 10.0};
   gsl_multiroot_function f = {&rosenbrock_f, n, &p};

   double x_init[2] = {-10.0, -5.0};
   gsl_vector *x = gsl_vector_alloc(n);

   gsl_vector_set(x, 0, x_init[0]);
   gsl_vector_set(x, 1, x_init[1]);

   T = gsl_multiroot_fsolver_hybrids;
   s = gsl_multiroot_fsolver_alloc(T, 2);

   gsl_multiroot_fsolver_set(s, &f, x);
   print_state(iter, s);

   do{
   iter++;
   status = gsl_multiroot_fsolver_iterate(s);
   print_state(iter, s);
   if(status){
   break;
   }
   status = gsl_multiroot_test_residual(s->f, 1e-7);
   }
   while(status == GSL_CONTINUE && iter < 1000);

   printf("status = %s\n", gsl_strerror(status));

   gsl_multiroot_fsolver_free(s);
   gsl_vector_free(x);
   return(0);
   }

//parse command line args
//create appropriate array
//read data into array
//generate unique failure times

*/
int count_lines(char *filename){

	int lines = 0;
	char ch;
	FILE *fp = fopen(filename, "r");

	while(!feof(fp)){
		ch = fgetc(fp);
		if(ch == '\n'){
			lines++;
		}
	}

	fclose(fp);
	return(lines);
}

struct Data{
	double *t;  //event time
	int *e;     //event indicator
	double *x1; //covariate 1
	int *x2;    //covariate 2
};

void read_data(char *f, struct Data x){

	FILE *fp = fopen(f, "r");

	double a;
	int b;
	double c;
	int d;

	int i = 0;
	while(fscanf(fp, "%lf %d %lf %d", &a, &b, &c, &d) != EOF){
		x.t[i] = a;
		x.e[i] = b;
		x.x1[i] = c;
		x.x2[i] = d;
		i++;
	}
	fclose(fp);
}

//sort all data by event time, in ascending order
//uses bubblesort, which is super inefficient
void bsort(struct Data data, int n){
	int i, newn, tempd;
	double templf;
	
	do{
		newn = 0;

		for(i=1; i < n; i++){
			if(data.t[i-1] > data.t[i]){

				//sort times in ascending order
				templf = data.t[i-1];
				data.t[i-1] = data.t[i];
				data.t[i] = templf;

				//sort event indicators by time
				tempd = data.e[i-1];
				data.e[i-1] = data.e[i];
				data.e[i] = tempd;

				//sort x1 by time
				templf = data.x1[i-1];
				data.x1[i-1] = data.x1[i];
				data.x1[i] = templf;

				//sort x2 by time
				tempd = data.x2[i-1];
				data.x2[i-1] = data.x2[i];
				data.x2[i] = tempd;

				newn = i;
			}
		}
		n = newn;
	} while( n != 0);
}

int main(int argc, char *argv[]){

	static char *in;
	static int n;

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

			//count number of observations
			n = count_lines(in);

			//create data structure 
			struct Data data;
			data.t = (double *)malloc(sizeof(double) * n);
			data.e = (int *)malloc(sizeof(int) * n);
			data.x1 = (double *)malloc(sizeof(double) * n);
			data.x2 = (int *)malloc(sizeof(int) * n);

			//read data, store in structure
			read_data(in, data);

			//sort by time
			bsort(data, n);

		}
	}
}
