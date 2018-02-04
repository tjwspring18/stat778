#include<stdlib.h>
int main(){
	/*float a[2619560];*/
	float *a = malloc(2619560 * sizeof(float));
	int i = 0;
	for(i; i < 2619560; i++){
		a[i] = 1.0;
	}

}
