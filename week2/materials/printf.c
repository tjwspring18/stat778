#include <stdio.h>
#include <string.h>

int main()
{

	printf("%.1f %.1f \n", 1.56, 1.25);
 
	printf("%x\n", 10);
	printf("%u\n", 10);
	printf("%c\n", 'A');
	printf("%f\n", 2.3);
	printf("%.1e\n", 0.012);
	printf("%d %%\n", 10);	

	printf("%4d\n",10);
	printf("%4s\n","hello");

	printf("%-4s","hello");

	char str[] = {'h','e','l','l','o','\0'};
	printf("%s\n", str);

	printf("%.2f & %.2f \\\\\n", 1.56, 1.25);
	printf("%d\n", strcmp("A","a"));
    
    double x;
    
    x = 0.1;
    
    printf("x=%f\n", x);
    
	return 0;
}
