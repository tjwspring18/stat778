/ . ensure included only once . / 
#ifndef __EUCLID_H__ 
#define __EUCLID_H__

/ . global variables (declared in euclid.c) . /
extern int x, y;

/ . compute gcd . /
int gcd ( int a, int b);

/ . compute g = gcd(a,b) and solve ax+by=g . /
int ext _euclid ( int a, int b);

#endif
