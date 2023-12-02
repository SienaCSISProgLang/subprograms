/*
    General purpose selection sort implementation, header file

    The selsort function is intended to be used interchangeably with
    the qsort, mergesort, and heapsort functions from the C standard library.

    Jim Teresco, Siena College, CSIS-340, Fall 2023
*/

#include <stdlib.h>

void selsort(void *base, size_t nmemb, size_t size, 
    int (*compar)(const void *, const void *));
