#include <stdio.h>

void xtest()
{
	const char str[] = "some text";

	if (str)
		printf("string is not null\n");
	else printf("string is null\n");

	///////////
	int x, y=0;
	for (x=1; x<=20; x++)
	{
		y = 0;
		if (x % 2 == 0)
			y += x/2;
		else if (x % 4 == 1)
			y += 2*((x+3)/4);
		else 
			y += (x+1)/2;
		printf("%d %d\n", x, y);	 
	}	

	for (x=1; x<=8; x++)
	{
		y = 0;
		if (x % 4 == 0)
			if (x % 2 == 0) 
				y = 2;
		else y = 1;
		printf("%d %d\n", x, y);
	}
}

void main()
{
	xtest();
}
