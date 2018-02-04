#include <stdio.h>

int main()
{
	int x=017; int y=12;
	
	if (x > y) printf("x > y\n");	

	int s=0xFFFF12; // short int 16 bits in size
	printf("%d\n", s);

	puts("hel""lo");

	//enum sz{S=0,L=3,XL};
	//printf("%d\n", XL);

	enum sz{S=0,L=-3,XL};
	printf("%d\n", XL);

	return 0;

}
