#include<stdio.h>
#include<stdlib.h>
#include "misc.h"

/*
   Counts number of lines in input file
   Used to know what size of arrays to create
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

/*
   Reads data from input file
   Stores in two arrays
   One for observation times (double)
   One for censor indicator (int)
 */
void read_data(char *filename, double A[], int B[]){
	FILE *fp = fopen(filename, "r");
	double a;
	int b;
	int i = 0;
	while(fscanf(fp, "%lf %d", &a, &b) != EOF){
		A[i] = a;
		B[i] = b;
		i++;
	}
	fclose(fp);
}

/* 
   Compares two values
   Is a smaller, larger, or same as b?
   Used in qsort function
 */
int cmpfunc(const void *a, const void *b){

	double x, y;

	x = *(double*)a;
	y = *(double*)b;

	if(x > y){
		return(1);
	}
	if(x == y){
		return(0);
	}
	if(x < y){
		return(-1);
	}

}

/* 
   Counts number of observed failures
 */
int count_failures(int B[], int n){
	int n_failures = 0;
	int i;
	for(i=0;i<n;i++){
		if(B[i] == 1){
			n_failures++;
		}
	}
	return(n_failures);
}

/* 
   Count unique values in unsorted array
   Used to know what size of array to create 
 */
int count_unique(double Af[], int nf){

	int i;
	int unique = 1;

	for(i=0; i < nf-1; i++){
		if(Af[i] == Af[i+1]){
			continue;
		} 
		else{
			unique++;
		}
	}
	return(unique);
}

void fill_unique_array(double Af[], double U[], int nfu){

	int i = 1;
	int j = 0;

	U[0] = Af[0]; //first element always unique

	while(i < nfu){
		if(Af[i] != U[j]){
			U[++j] = Af[i];
		}
		i++;
	}
}

//count number of failures at each unique failure time
void count_failed(double A[], int B[], int D[], 
		double U[], int n, int nfu){
	int i, j, d;
	for(i=0; i < nfu; i++){
		d = 0;
		for(j=0; j < n; j++){
			if(U[i] == A[j]){
				d = d + B[j];
			}
		}
		D[i] = d;
	}
}

void calculate_risk_set(double A[], int R[], double U[], 
		int n, int nfu){

	int i, j, r;
	for(i=0; i<nfu; i++){
		r = n;
		for(j=0; j<n; j++){
			if(U[i] > A[j]){
				r--;
			} 
		}
		R[i] = r;;
	}
}

/* 
   Write data to output file
 */
void write_data(char *filename, 
		double U[], 
		int R[], 
		int D[],
		double S[],
		double V[],
		double CIL[],
		double CIU[],
		int n){

	FILE *fp = fopen(filename, "w");
	int i = 0;

	//fprintf(fp, "%s,%s,%s\n", "lower", "s(t)", "upper");
	fprintf(fp, "%s,%s,%s,%s,%s,%s,%s,\n",
			"time",
			"risk_set",
			"failures",
			"survival",
			"se",
			"lower95", 
			"upper95");

	for(i; i < n; i++){
		fprintf(fp, "%lf,%d,%d,%lf,%lf,%lf,%lf\n", 
				U[i],
				R[i],
				D[i],
				S[i],
				V[i],
				CIL[i],
				CIU[i]);
	}
	fclose(fp);
}
