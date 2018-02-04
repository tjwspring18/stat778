// routines for calculating poission distributions

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
    
    prob  = pow(mu, x)*exp(-mu)/factorial1(x);
    
    return prob;
}

int factorial1(int n)
{
	int i, j=1;
	for (i=1; i<=n; i++)
		j *= i;
	return j;
}