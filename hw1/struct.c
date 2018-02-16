#include<stdio.h>
typedef struct my_data{
	int a[2];
	int b[2];
} my_data;

void main(){
	//works
	my_data foo;
	foo.a[0] = 1;
	printf("%d\n", foo.a[0]);
	/*
	struct my_data data[2];

	static int data_cmp(const void *a, const void *b){
		const struct my_data *da = a, *db = b;
		return da->a < db->a ? -1 : da->a > db->a;
	}

	*/
}
