int count_lines(char *filename);
void read_data(char *filename, double A[], int B[]);
int cmpfunc(const void *a, const void *b);
int count_failures(int B[], int n);
int count_unique(double Af[], int nf);
void fill_unique_array(double A[], double U[], int n);
void count_failed(double A[], int B[], int D[], 
		double U[], int n, int nfu);
void calculate_risk_set(double A[], int R[], double U[], 
		int n, int nfu);
void write_data(char *filename, 
		double U[], 
		int R[], 
		int D[],
		double S[],
		double V[],
		double CIL[],
		double CIU[],
		int n);
