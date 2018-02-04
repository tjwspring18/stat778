/*
   Tom Wallace <twalla11@masonlive.gmu.edu>

   STAT 778
   Spring 2018
   HW #1 submission: Kaplan-Meier estimator
*/

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include "misc.h"

int main(int argc, char *argv[]){
	
	//Print information about program
	printf("Running Kaplan-Meier estimator program\n");
	printf("Compiled %s %s \n", __DATE__, __TIME__);
	
	//variable declarations
	extern char *optarg;
	int iflag = 0, oflag = 0; 
	char *iname, *oname; 
	int c, i, j, q = 0;
	char *abort_message = "Program aborted before successful completion";
	FILE *fpi;
	int n;
	int nu;

	//parse command line arguments
	//required: -i input_file -o output_file
	while (( c = getopt(argc, argv, "i:o:") ) != -1){
		switch (c) {
			case 'i':
				iflag = 1;
				iname = optarg;
				break;
			case 'o':
				oflag = 1;
				oname = optarg;
				break;
		}
	}

	//abort if input and output files not received
	if( (iflag = 0) || (oflag == 0)){

		printf("Error: did not receive arguments -i input_file -o output_file\n");
		printf("%s\n", abort_message);
		exit(1);

	} else{
		//otherwise continue
		printf("Using input file: %s\n", iname);

		//check that input file exists
		//abort if it does not exist
		if( access(iname, F_OK) == -1 ){
			printf("Input file not found\n");
			printf("%s\n", abort_message);
			exit(1);
		} else{
			//count lines in input file
			n = count_lines(iname);

			//create 2 arrays of that length
			double A[n];
			int B[n];

			//read data from file into those arrays
			read_data(iname, A, B);
			printf("Read %d observations from input file\n", n);
			printf("Computing Kaplan-Meier estimates\n");
			
			//sort A as new array As
			double As[n];
			for(i; i < n; i++){
				As[i] = A[i];
			}
			qsort(As, 200, sizeof(double), cmpfunc);

			//count number of unique failure times in As
			nu = count_unique(As, n);

			//create new array of that length
			double U[nu];

			//copy unique failure times from A to U
			fill_unique_array(As, U, n);

			//create array C of length nu 
			//calculate number censored at each unique time
			int C[nu];
			for(i=0; i < nu; i++){
				q = 0;
				for(j=0; j < n; j++){
					if(U[i] == A[j]){
						q = q + (1-B[j]);
					}
				}
				C[i] = q;
			}
			

			//create array D of length nu
			//calculate number died at each unique time
			exit(0);
		}
	}
}
