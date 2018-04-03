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
	double *x0; // = 1 (intercept term)
	double *x1; // covariate 1
	double *x2; // covariate 2
	size_t size;
};
struct Data* makeData(size_t sz);
void deleteData(struct Data *data);

int count_lines(char *f); 
void read_data(char *f, struct Data *data);
void bsort(struct Data *data, int n);

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

			/*
			for(i=0; i<n; i++){
				printf("%lf %d %lf %lf %lf\n",
						data->t[i],
						data->e[i],
						data->x0[i],
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
	data->x0 = malloc(sz * sizeof(double));
	data->x1 = malloc(sz * sizeof(double));
	data->x2 = malloc(sz * sizeof(double));
	data->size = sz;
	return(data);
}

// delete instance of Data
void deleteData(struct Data *data){
	free(data->t);
	free(data->e);
	free(data->x0);
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
		data->x0[i] = 1;
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
