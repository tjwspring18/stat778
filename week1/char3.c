#include <stdio.h>

/* count lines
in input
 */
int main()
{
    int c, nc;     /* one
                      two
             three  */
    nc = 0;
    while ((c=getchar())!=EOF)
        if (c == '\n')
            ++nc;;
    
    printf("%d\n", nc);
    return 0;
    
    
}
