#include <stdio.h>

int main()
{
	int x, y;

	for (x=1; x<=10; x++)
	{
		y = 0;
		if (x % 2)
			y+= x/2;
		printf("%d %d\n", x, y);
	}

	return 0;
}	
