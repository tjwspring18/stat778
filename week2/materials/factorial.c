#include <stdio.h>

int factorial1(int n)
{
	int i, j=1;
	for (i=1; i<=n; i++)
		j *= i;
	return j;
}

int factorial2(int n)
{
	int i=1, j=1;
	
    printf("Factorial 2:\n");
	while(i<=n)
	{
		j *= i;
		i++;
        printf("i=%d, j=%d\n", i, j);
	}
	return j;
}

// P(X=x), where X follows a poisson distribution
// with mean mu
double dpois(int x, double mu)
{
    double prob;
    
    prob  = pow(mu, x)*exp(-x)/factorial1(x);
    
    return prob;
}

int main()
{
	int n=5, j=6;

	printf("%d %d\n", factorial1(n), factorial2(n));
    
    factorial2(n);
    
    printf("factorial(%d)=%d\n",j,factorial1(j));
    
	return 0;
}
