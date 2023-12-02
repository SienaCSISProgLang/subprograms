/*
  A function to compute the greatest common denominator of two integers.

  Implementation uses the Euclidean Algorithm.

  Jim Teresco

  Fall 2009, Mount Holyoke College
  Fall 2010, Siena College
  Fall 2012, The College of Saint Rose

*/

/* gcd function: takes two integer parameters, uses Euclidean
   Algorithm to compute and return greatest common denominator.

   Requirement: a and b are nonnegative.
*/
int gcd(int a, int b) {

  if (a == 0) return b;
  if (a > b) return gcd(b, a);
  return gcd(b%a, a);
}
