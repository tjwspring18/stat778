#include <stdio.h>

// print Fahrenheit-Celsius table for fahr=0,20,...,300

int main()
{
    int fahr, Celsius;
    int lower, upper, step;
    
    lower = 0; // lower limit of temperature table
    upper = 300; // upper limit
    step = 20; // step size
    
    fahr = lower;
    while (fahr <= upper) {
        Celsius = 5*(fahr-32)/9;
        printf("%d\t%d\n", fahr, Celsius);
        fahr += step;
    }
    return 0;
}
