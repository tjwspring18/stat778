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

int main(int argc, char *argv[]){
	
	//Print information about program
	printf("Running Kaplan-Meier estimator program\n");
	printf("Compiled %s %s \n", __DATE__, __TIME__);
	
	//variable declarations
	extern char *optarg;
	int iflag = 0, oflag = 0; 
	char *iname, *oname; 
	int c, i, j, q = 0;
	double s;
	char *abort_message = "Program aborted before successful completion";
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
			double *A = malloc(n * sizeof(double));
			int *B = malloc(n * sizeof(int));

			//read data from file into those arrays
			read_data(iname, A, B);
			printf("Read %d observations from input file\n", n);
			printf("Computing point estimates\n");

			//sort A into new array As
			double *As = malloc(n * sizeof(double));
			for(i=0; i < n; i++){
				As[i] = A[i];
			}
			qsort(As, 200, sizeof(double), cmpfunc);

			//count number of unique failure times in As
			nu = count_unique(As, n);

			//create new array U of that length
			//copy unique failure times from As to U
			double *U = malloc(nu * sizeof(double));
			fill_unique_array(As, U, n);

			//dont need As anymore
			free(As);

			//create array C of length nu 
			//calculate number censored at each unique time
			int *C = malloc(nu * sizeof(int));
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
			int *D = malloc(nu * sizeof(int));
			for(i=0; i < nu; i++){
				q = 0;
				for(j=0; j < n; j++){
					if(U[i] == A[j]){
						q = q + B[j];
					}
				}
				D[i] = q;
			}
			
			//dont need A or B anymore
			free(A);
			free(B);

			//create array R of length nu
			//calculate risk set at each unique time
			int *R = malloc(nu * sizeof(int));
			for(i = 0; i < nu; i++){
				q = 0;
				for(j = 0; j < i; j++){	
					q = q + D[j] + C[j];
				}
				R[i] = n-q;
			}

			//marginal survival for each time
			double *Sm = malloc(nu * sizeof(double));
			for(i = 0; i < nu; i++){
				Sm[i] = ((double)R[i] - (double)D[i]) / (double)R[i];
			}


			//cumulative survival for each time 
			double *S = malloc(nu * sizeof(double));
			for(i = 0; i < nu; i++){
				s = 1;
				for(j = 0; j < i; j++){	
					s = s * Sm[j];
				}
				S[i] = s;
			}

			printf("Computing variance\n");

			//stddev for each S(t)
			//something is going wrong here 
			//maybe has to do with D[i] being 0 for censored observations
			//unique times should be only be for unique DEATH times
			double *V = malloc(nu * sizeof(double));
			for(i = 0; i < nu; i++){
				s = 0;
				for(j = 0; j < i; j++){	
					s = s + ((double)D[i] / ( (double)R[i] * ((double)R[i] - (double)D[i]) ));
				}
				V[i] = sqrt(s * pow(S[i], 2));
			}

			//95 percent CI, lower 
			double *CIL = malloc(nu * sizeof(double));
			for(i = 0; i < nu; i++){
				CIL[i] = S[i] - V[i]*1.96;
			}

			//95 percent CI, upper
			double *CIU = malloc(nu * sizeof(double));
			for(i = 0; i < nu; i++){
				CIU[i] = S[i] + V[i]*1.96;
			}

			//write results
			printf("Writing to output file: %s\n", oname);
			write_data(oname, CIL, S, CIU, nu);

			exit(0);
		
		}
	}
}
