#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#define ABORT_MESSAGE "Program aborting"
#define BAD_ARGS "Wrong number of arguments passed"
#define FNF "File not found"


struct Data{
	double *t;  //event time
	int *e;     //event indicator
	double *x1; //covariate 1
	double *x2;    //covariate 2
};
int count_lines(char *filename); //count number of observations
void read_data(char *f, struct Data x); //read input file into Data structure
void bsort(struct Data data, int n); //sort of Data by time


int main(int argc, char *argv[]){

	static char *in;
	static int n;
	int i, j;
	double b1 = 0.5;
	double b2 = -0.5;
	double l, q, r=0;

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
			printf("Read %d observations\n", n);

			//create data structure 
			struct Data data;
			data.t = (double *)malloc(sizeof(double) * n);
			data.e = (int *)malloc(sizeof(int) * n);
			data.x1 = (double *)malloc(sizeof(double) * n);
			data.x2 = (double *)malloc(sizeof(double) * n);

			//read data, store in structure
			read_data(in, data);

			//sort by time
			bsort(data, n);
			
			//calculate partial log likelihood with
			//beta1=0.5 and beta2=-0.5
			printf("Calculating log likelihood\n");
			for(i=0; i<n; i++){

				if(data.e[i] == 1){

					q = b1*data.x1[i] + b2*data.x2[i];

					for(j=i; j<n; j++){
						r += exp(b1*data.x1[j] + b2*data.x2[j]);
					}

					l += q - log(r);
				}
			}
			printf("Log likelihood with beta1=0.5 and beta2=-0.5: %lf\n", l);
		}
	}
}

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

void read_data(char *f, struct Data x){

	FILE *fp = fopen(f, "r");

	double a;
	int b;
	double c;
	double d;

	int i = 0;
	while(fscanf(fp, "%lf %d %lf %lf", &a, &b, &c, &d) != EOF){
		x.t[i] = a;
		x.e[i] = b;
		x.x1[i] = c;
		x.x2[i] = d;
		i++;
	}
	fclose(fp);
}

//this uses bubblesort, which is computationaly inefficient (but conceptually
//simple and hence easy for me to code)
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
