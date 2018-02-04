#include <stdio.h> // basic I/O facilities
#define msg2 "hello, world"

// The main() function
int main()
{
	const char msg[] = "hello, world";

	// write message to console
	printf("hello, world\n");
	
	puts("hello, world\n");

	printf("%s\n", msg);

	puts(msg);

	puts(msg2);
	
	return 0; // exit (0 = success)
}
