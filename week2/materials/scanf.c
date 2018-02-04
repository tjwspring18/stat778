#include <stdio.h>

int main()
{
	int x;	
	double y;
	char str[100], fname[100], lname[100];

	scanf("%d", &x);
	printf("%d\n", x);

	scanf("%lf", &y);
	printf("%f\n", y);

	scanf("%s",str);
	printf("%s\n", str);

	scanf("%20s %20s", fname, lname);
	printf("%s %s\n", fname, lname);

	return 0;
}
