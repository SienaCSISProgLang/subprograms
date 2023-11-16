/*
  A slightly more interesting but still very simple C program, which
  accepts two numbers as command-line parameters or entered by the
  user when prompted, and computes and prints out their greatest
  common denominator.

  Jim Teresco

  Fall 2010, Siena College
  Fall 2012, The College of Saint Rose

  $Id$
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>

#include "gcd.h"

int main(int argc, char *argv[]) {
  int a, b, retval;

  if (argc < 3) {
    printf("Please enter the two numbers whose gcd you wish to compute: ");
    retval = scanf("%d%d", &a, &b);
    if (retval != 2) {
      fprintf(stderr, "%s: Could not parse two integer values from input.\n",
	      argv[0]);
      exit(1);
    }
  }
  else {
    a = strtol(argv[1], NULL, 10);
    if (errno) {
      fprintf(stderr, "%s: Could not convert %s to an integer value\n", 
	      argv[0], argv[1]);
      exit(1);
    }
    b = strtol(argv[2], NULL, 10);
    if (errno) {
      fprintf(stderr, "%s: Could not convert %s to an integer value\n", 
	      argv[0], argv[2]);
      exit(1);
    }
  }

  if ((a < 0) || (b < 0)) {
    fprintf(stderr, "%s: Both numbers must be non-negative\n", argv[0]);
    exit(1);
  }

  printf("GCD of %d and %d is %d\n", a, b, gcd(a,b));

  return 0;
}
