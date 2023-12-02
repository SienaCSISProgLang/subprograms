/*
  CS 322/324 Example Code
  Mount Holyoke College
  CS 220, Siena College
  CSC 433, The College of Saint Rose
  Jim Teresco

  Header file for ratio structure and related functions

*/

#ifndef __RATIO_H
#define __RATIO_H

typedef struct ratio {
  int numerator;
  int denominator;
} ratio;

extern ratio *create_ratio(int numerator, int denominator);
extern ratio *add_ratios(ratio *rat1, ratio *rat2);
extern void reduce_ratio(ratio *rat);
extern double ratio_to_double(ratio *rat);
extern void print_ratio(ratio *rat);
extern int ratio_numerator(ratio *rat);
extern int ratio_denominator(ratio *rat);

#endif
