/*
  Demonstation of Jensen's Device using C preprocessor

  Based on example from
  http://ww2.valdosta.edu/~dboyd/CS4330/Lectures/Jensen.pdf

  Updated by Jim Teresco, The College of Saint Rose, Fall 2014, CSC 433
 */

#include <stdio.h>
#include <math.h>

/* macro to implement Jensen's device */

#define SUM(Expr,Term,Lbound,Ubound,Inc,Total) \
  for(Total = 0, Term = Lbound; Term <= Ubound; Term += Inc)\
    Total += (Expr)

int main() {
  int X;
  float a[]={1,2,3,4,5,6,7}, b[]={1,-2,3,-4,5,-6,7};
  float T;

  SUM(X,X,1,100,2,T);
  printf("Sum odds 1 thru 99: %g\n",T);

  SUM(1.0/a[X],X,0,6,1,T);
  printf("Sum of reciprocals: %g\n",T);

  SUM(a[X]*b[X],X,0,6,1,T);
  printf("Dot product: %g\n",T);

  return 0;
}

