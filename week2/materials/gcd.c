#include <stdio.h>

#define min(a,b) ((a) < (b) ? (a):(b))
#define max(a,b) ((a) < (b) ? (b):(a))

/* common divisor of a and b */
int gcd(int a, int b)
{
	int i, ret=1, minval=min(a,b);

	for (i=2; i<=minval; i++)
	{
		if (a % i) /* i not divisor of a */
			continue;
		if (b % i == 0) /* i is divisor of both a and b */
			ret = i;
	}

	return ret;
}

int gcd2(int a , int b)
{
	int a1, b1;

	a1 = max(a,b);
	b1 = min(a,b);

	while (b1)
	{
		int temp=b1;
		b1 = a1 % b1;
		a1  = temp;
	} 
	return a1;
}

int main()
{
	int a, b;

	for (a=1; a<=5; a++)
		for (b=1; b<=5; b++)
			printf("%d %d %d %d\n", a, b, gcd(a,b), gcd2(a,b));
	return 0;
}
