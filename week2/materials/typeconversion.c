#include <stdio.h>

int main()
{
    int i, j;
    double f;
    
    f = 1.75;
    i = (int)f;
    
    printf("%f %d\n", f, i);
    
    i = 3; j = 4;
    
    printf("%d %d %f\n", i, j, 1.0*i/j);
    
}

