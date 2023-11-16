/*
  Variable argument list length example in C

  Jim Teresco
  The College of Saint Rose
  Programming Languages, CSC 433
  Fall 2012

*/

#include <stdio.h>
// we will need this header file if we plan to use variable argument lists
#include <stdarg.h>

// variable length parameter list of ints to add up -- must end in a 
// negative number as a sentinel
int sumnumbers(int first, ...) {

  va_list argp;  // argument pointer

  va_start(argp, first); // initialize to point to "first"

  // now loop over the arguments and accumulate the sum
  int sum = 0;
  int number;
  for (number = first; number >= 0; number = va_arg(argp, int)) {
    sum += number;
  }
  va_end(argp);

  return sum;
  
}

int main() {

  printf("First: %d\n", sumnumbers(1, 2, 3, -1));
  printf("Second: %d\n", sumnumbers(82, 12, 21, 211, 98, 7, 12, 87, -23));
  printf("Third: %d\n", sumnumbers(-100));
  // what if we don't have a sentinel??
  printf("Fourth: %d\n", sumnumbers(1, 2, 3, 4, 5, 6, 7, 8, 9, 10));
  // what if we give it the wrong parameters?
  printf("Fifth: %d\n", sumnumbers(1, 8.2, "hi there!", 87L, -1));
  return 0;
}
