/*
   Tom Wallace <twalla11@masonlive.gmu.edu>

   STAT 778
   Spring 2018
   HW #1 submission: Kaplan-Meier estimator
*/
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include "misc.h"

#define ABORT_MESSAGE "Program aborted before successful completion"

int main(int argc, char *argv[]){
	
	//Print information about program
	printf("Running Kaplan-Meier estimator program\n");
	printf("Compiled %s %s \n", __DATE__, __TIME__);
	
	//variable declarations
	extern char *optarg;
	int iflag = 0, oflag = 0; 
	char *iname, *oname; 
	int c, i, j, q = 0;
	int n, nf, nfu;
	double s;

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
		printf("%s\n", ABORT_MESSAGE);
		exit(1);

	} else{
		printf("Using input file: %s\n", iname);

		//check that input file exists
		//abort if it does not exist
		if( access(iname, F_OK) == -1 ){
			printf("Input file not found\n");
			printf("%s\n", ABORT_MESSAGE);
			exit(1);
		} else{
			//count lines in input file
			n = count_lines(iname);

			//create 2 arrays of that length
			double *A = malloc(n * sizeof(double));
			int *B = malloc(n * sizeof(int));

			//read data from file into those arrays
			read_data(iname, A, B);
			printf("Read %d observations from input file\n", n);
			printf("Computing point estimates\n");

			//count number of failure time observations
			//not necessarily n of UNIQUE failure time observations
			nf = count_failures(B, n);

			//create new array Af of length nf
			//fill it with failure time observations
			double *Af = malloc(nf * sizeof(double));

			for(i=0; i < n; i++){
				if(B[i] == 1){
					Af[q] = A[i];
					q++;
				}
			}

			//sort Af 
			//count number of unique failure times in Af
			qsort(Af, nf, sizeof(double), cmpfunc);
			nfu = count_unique(Af, nf);

			//create new array U of that length
			//copy unique failure times from Af to U
			//double check that this works by modifying input data
			double *U = malloc(nfu * sizeof(double));
			fill_unique_array(Af, U, nfu);

			//dont need Af anymore
			free(Af);

			//create array D of length nfu
			//calculate number failed at each unique failure time
			int *D = malloc(nfu * sizeof(int));
			count_failed(A, B, D, U, n, nfu);

			//create array R of length nfu 
			//calculate risk set at each unique failure time
			//store in R
			int *R = malloc(nfu * sizeof(int));
			calculate_risk_set(A, R, U, n, nfu);

			//dont need A or B anymore
			free(A);
			free(B);

			//marginal survival for each time
			double *Sm = malloc(nfu * sizeof(double));
			for(i = 0; i < nfu; i++){
				Sm[i] = ((double)R[i] - (double)D[i]) / (double)R[i];
			}


			//cumulative survival for each time 
			double *S = malloc(nfu * sizeof(double));
			for(i = 0; i < nfu; i++){
				s = 1.0;
				for(j = 0; j < i; j++){	
					s *= Sm[j];
				}
				S[i] = s;
			}

			//get rid of Sm
			free(Sm);

			printf("Computing variance\n");

			// Variance, Greenwood's formula
			//
			// Definitions
			//
			// V[i] = variance of KM estimator (what I am computing here)
			// nfu  = number of unique failure times 
			// D[i] = number of failures at time t_i
			// R[i] = risk set at time t_i
			// S[i] = KM estimator at time t_i 
			double *V = malloc(nfu * sizeof(double));
			for(i = 0; i < nfu; i++){
				s = 0;
				for(j = 0; j < i; j++){	
					s = s + (((double)D[i]) / ( (double)R[i] * ((double)R[i] - (double)D[i]) ));
				}
				V[i] = s * pow(S[i], 2);
			}

			//standard error of KM estimator
			//square root of variance calculated above
			double *SE = malloc(nfu * sizeof(double));
			for(i=0; i < nfu; i++){
				SE[i] = sqrt(V[i]);
			}

			//95 percent CI, lower 
			double *CIL = malloc(nfu * sizeof(double));
			for(i = 0; i < nfu; i++){
				CIL[i] = pow(S[i], exp( (-1.96*SE[i]) / (S[i]*log(S[i])) ));
			}

			//95 percent CI, upper
			double *CIU = malloc(nfu * sizeof(double));
			for(i = 0; i < nfu; i++){
				CIU[i] = pow(S[i], exp( (1.96*SE[i]) / (S[i]*log(S[i])) ));
			}

			//write results
			printf("Writing to output file: %s\n", oname);
			write_data(oname, U, R, D, S, SE, CIL, CIU, nfu);

			//clean up
			free(CIU);
			free(CIL);
			free(S);
			free(V);
			free(U);
			free(D);
			free(R);
			free(SE);

			exit(0);
		
		}
	}
}
