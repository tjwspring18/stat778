#include <stdio.h>

int main()
{
	char c;
	int x = 0;

	do {
		puts("Keep going? (y/n) ");
		c = getchar();
		//printf("%c\n", c);
		if (c != 'y' && c != '\n')
        //if (c != 'y')
			break;
		printf("continue\n");
	} while (x >= 0);

	return 0;
}
