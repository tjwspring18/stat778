#include <stdio.h>

// print Fahrenheit-Celsius table for fahr=0,20,...,300

int main()
{
    double fahr, Celsius;
    
    for (fahr=0; fahr<=300; fahr = fahr+20)
        printf("%.0f %6.1f\n", fahr, (5.0/9.0)*(fahr-32));
    return 0;
}
