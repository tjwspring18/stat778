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
	int unique = 0;

	for(i=0; i < nf; i++){
		if(Af[i] == Af[i+1]){
			continue;
		} 
		else{
			unique++;
		}
	}
	return(unique);
}

void fill_unique_array(double A[], double U[], int n){

	int i = 1;
	int j = 0;

	U[0] = A[0]; //first element always unique

	while(i < n){
		if(A[i] != U[j]){
			U[++j] = A[i];
		}
		i++;
	}
}

void count_censored(double A[], int B[], double U[], n, nf)
/* 
   Write data to output file
*/
void write_data(char *filename, double *CIL, double *S, double *CIU, int n){
	
	FILE *fp = fopen(filename, "w");
	int i = 0;
	
	fprintf(fp, "%s,%s,%s\n", "lower", "s(t)", "upper");
	
	for(i; i < n; i++){
		fprintf(fp, "%lf,%lf,%lf\n", CIL[i], S[i], CIU[i]);
	}
	fclose(fp);
}
