/*
  Example of C pointer parameters acting like reference parameters

  Jim Teresco
  The College of Saint Rose
  Programming Languages, CSC 433
  Fall 2012
*/

#include <stdio.h>

// both x and y are passed by value, but y is a pointer to a variable
// in the caller
void f(int x, int *y) {

  x = 23;
  *y = 17;
}


int main(int argc, char *argv[]) {

  int a, b;
  a = 1;
  b = 2;

  printf("Before: a=%d, b=%d\n", a, b);

  f(a, &b);

  printf("After: a=%d, b=%d\n", a, b);

  // this won't be allowed, since "2" can't be changed
  //f(1, &2);

  return 0;
}
