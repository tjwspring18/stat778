int count_lines(char *filename);
void read_data(char *filename, double A[], int B[]);
int cmpfunc(const void *a, const void *b);
int count_failures(int B[], int n);
int count_unique(double Af[], int nf);
void fill_unique_array(double A[], double U[], int n);
void write_data(char *filename, double *CIL, double *S, double *CIU, int n);
