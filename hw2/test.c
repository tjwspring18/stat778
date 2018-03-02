#include<stdlib.h>
#include<stdio.h>

int main(void){

	int a[] = {3, 1, 2};
	int b[] = {3, 1, 2};
	int temp, i, newn;
	int n;

	n = sizeof(a) / sizeof(int);

	do{
		newn = 0;
		for(i=1; i < n; i++){
			if(a[i-1] > a[i]){
				temp = a[i-1];
				a[i-1] = a[i];
				a[i] = temp;
				temp = b[i-1];
				b[i-1] = b[i];
				b[i] = temp;
				newn = i;
			}
		}
		n = newn;
	} while( n != 0);

	for(i=0; i < 3; i++){
		printf("%d %d\n", a[i], b[i]);
	}

}
