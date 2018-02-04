#include <stdio.h>

int main(int argc, char **argv)
{
	char  x[100];
    double y, z;

	FILE *fp, *fs;

	if ((fp=fopen(argv[1],"r"))==NULL)
		printf("Fail to read the input file %s!\n", argv[1]);
	
	fscanf(fp, "%s", x);
    fscanf(fp, "%lf%lf", &y, &z);
	fclose(fp);

	if ((fs=fopen(argv[2],"w"))==NULL)
                printf("Fail to open output file %s!\n", argv[2]);
	fprintf(fs, "The content of %s: %s\n", argv[1], x);
    fprintf(fs, "%f %f\n", y, z);
	fclose(fs);
	
	return 0;
}
