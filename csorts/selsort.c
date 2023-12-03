/*
    General purpose selection sort implementation

    The selsort function is intended to be used interchangeably with
    the qsort, mergesort, and heapsort functions from the C standard library.

    Jim Teresco, Siena College, CSIS-340, Fall 2023
    Aided by Copilot
*/

#include <string.h>
#include "selsort.h"

void selsort(void *base, size_t nmemb, size_t size,
         int (*compar)(const void *, const void *)) {

  // we'll need a temporary buffer to hold the smallest element
  void *temp = malloc(size);

  // outer loop: we'll be filling in the sorted portion of the array
  for (int i = 0; i < nmemb; i++) {

    // find the smallest element in the unsorted portion of the array
    int smallest = i;
    for (int j = i+1; j < nmemb; j++) {
      if (compar(base + j*size, base + smallest*size) < 0) {
        smallest = j;
      }
    }

    // swap the smallest element with the first element of the unsorted
    // portion of the array
    memcpy(temp, base + i*size, size);
    memcpy(base + i*size, base + smallest*size, size);
    memcpy(base + smallest*size, temp, size);
  }

  // free the temporary buffer
  free(temp);
}

// for comparison, a simpler selection sort implementation that works only on
// arrays of int and only sorts in ascending order
void selsort_int(int *base, size_t nmemb) {

  // outer loop: we'll be filling in the sorted portion of the array
  for (int i = 0; i < nmemb; i++) {

    // find the smallest element in the unsorted portion of the array
    int smallest = i;
    for (int j = i+1; j < nmemb; j++) {
      if (base[j] < base[smallest]) {
        smallest = j;
      }
    }

    // swap the smallest element with the first element of the unsorted
    // portion of the array
    int temp = base[i];
    base[i] = base[smallest];
    base[smallest] = temp;
  }
}
