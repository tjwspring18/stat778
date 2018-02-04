#include <stdio.h>

// copy input to output
int main()
{
    int c;
    
    /*
    c = getchar();
    while (c != EOF) {
        putchar(c);
        c=getchar();
    }*/
    
    while ((c=getchar())!=EOF)
        putchar(c);
    
    //while ((c=getchar())!='c') {
    //    putchar(c);
    //}
    return 0;
    
    
}
