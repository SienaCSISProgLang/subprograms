/*
  Demonstrations of general-purpose sorting in C
  
  Based on: sxamples of usage of the qsort function from the C standard library.
  Jim Teresco, CSC 381, The College of Saint Rose, Fall 2013

  Initial implementation:
  Sat Nov  9 13:08:41 EST 2013

  Expanded to other sorting for
  CSIS-340, Siena College, Fall 2023
*/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include "ratio.h"
#include "selsort.h"

/* types of data we can sort */
enum datatypes { NUMERIC, STRING, RATIO };
/* special ways we might sort them */
enum criteria { DEFAULT, LENGTH, CASEINSENSITIVE, NUMERATOR };

/* should allocate dynamically, but this makes life easier for this example */
#define MAX_ELEMENTS 10000
#define MAX_WORDLEN 128

/* it's always nice to have a usage function to call when you find out
   that something isn't right on the command line */
static void usage(char *program) {
  printf("Usage: %s [--filename filename] [--datatype numeric|string|ratio] [--reverse]\n", program);
  printf("[--criteria length|case-insensitive|numerator] [--algorithm selection|merge|quick|heap]\n");
  printf("[--verbose] [--help]\n");
}

/***********************************************************************/
/*  COMPARATOR FUNCTIONS                                               */
/***********************************************************************/

/* compare two int values, standard order */
int intcmp(const void *i1, const void *i2) {

  return *(int *)i1 - *(int *)i2;
}


/* compare two int values, reverse order */
int intcmprev(const void *i1, const void *i2) {

  return *(int *)i2 - *(int *)i1;
}

/* compare two int values, number of digits only */
int intcmpdigits(const void *i1, const void *i2) {

  int int1 = *(int *)i1;
  int int2 = *(int *)i2;

  // count up the digits
  int digits1 = 0;
  int digits2 = 0;

  while (int1) {
    digits1++;
    int1 /= 10;
  }

  while (int2) {
    digits2++;
    int2 /= 10;
  }

  return digits1 - digits2;
}

/* compare two int values, reverse number of digits */
int intcmpdigitsrev(const void *i1, const void *i2) {

  return -intcmpdigits(i2, i1);
}

/* wrapper for standard string compare */
int stringcmp(const void *s1, const void *s2) {

  char *str1 = *(char **)s1;
  char *str2 = *(char **)s2;

  return strcmp(str1, str2);

}


/* wrapper for standard string compare, case insensitive */
int stringcmpnocase(const void *s1, const void *s2) {

  char *str1 = *(char **)s1;
  char *str2 = *(char **)s2;

  return strcasecmp(str1, str2);

}

/* reverse of standard string compare */
int stringcmprev(const void *s1, const void *s2) {

  char *str1 = *(char **)s1;
  char *str2 = *(char **)s2;

  return strcmp(str2, str1);

}

/* and reverse no case */
int stringcmpnocaserev(const void *s1, const void *s2) {

  char *str1 = *(char **)s1;
  char *str2 = *(char **)s2;

  return strcasecmp(str2, str1);

}

/*  string compare by length */
int stringcmplength(const void *s1, const void *s2) {

  char *str1 = *(char **)s1;
  char *str2 = *(char **)s2;

  return strlen(str1) - strlen(str2);

}

/*  string compare by length */
int stringcmplengthrev(const void *s1, const void *s2) {

  char *str1 = *(char **)s1;
  char *str2 = *(char **)s2;

  return strlen(str2) - strlen(str1);

}

/* standard ratio comparison */
int ratiocmp(const void *rat1, const void *rat2) {

  ratio *r1 = *(ratio **)rat1;
  ratio *r2 = *(ratio **)rat2;

  double r1val = ratio_to_double(r1);
  double r2val = ratio_to_double(r2);

  if (r1val < r2val) return -1;
  if (r1val == r2val) return 0;
  return 1;
}


/* reverse ratio comparison */
int ratiocmprev(const void *rat1, const void *rat2) {

  ratio *r1 = *(ratio **)rat1;
  ratio *r2 = *(ratio **)rat2;

  double r1val = ratio_to_double(r1);
  double r2val = ratio_to_double(r2);

  if (r1val < r2val) return 1;
  if (r1val == r2val) return 0;
  return -1;
}


/* numerator-only ratio comparison */
int ratiocmpnum(const void *rat1, const void *rat2) {

  ratio *r1 = *(ratio **)rat1;
  ratio *r2 = *(ratio **)rat2;

  return ratio_numerator(r1) - ratio_numerator(r2);
}

/* numerator-only reverse ratio comparison */
int ratiocmpnumrev(const void *rat1, const void *rat2) {

  ratio *r1 = *(ratio **)rat1;
  ratio *r2 = *(ratio **)rat2;

  return ratio_numerator(r2) - ratio_numerator(r1);
}

/**********************************************************************/

int main(int argc, char *argv[]) {

  /* a file name -- empty to start, and we'll use stdin if none is specified */
  char filename[FILENAME_MAX] = "";

  /* the file pointer we'll read from */
  FILE *infp;

  /* other params we'll set */
  enum datatypes datatype = NUMERIC;
  enum criteria criterion = DEFAULT;

  /* variable needed for getopt_long */
  int ch;

  /* define our "long" equivalents for the short options */

  struct option longopts[] = {
    { "filename", required_argument, NULL, 'f' },
    { "datatype", required_argument, NULL, 'd' },
    { "reverse", no_argument, NULL, 'r' },
    { "criteria", required_argument, NULL, 'c' },
    { "algorithm", required_argument, NULL, 'a'},
    { "help", no_argument, NULL, 'h' },
    { "verbose", no_argument, NULL, 'v' },
    { 0, 0, 0, 0 }
  };

  int option_index = 0;

  /* was verbose specified?  Default of 0 will be overridden by -v */
  int verbose = 0;

  /* reverse sort? */
  int reverse = 0;

  /* storage for the actual arrays we might use */
  int numbers[MAX_ELEMENTS];
  char *words[MAX_ELEMENTS];
  ratio *fractions[MAX_ELEMENTS];

  char word[MAX_WORDLEN];
  int item;
  int num, den;

  /* pointer to the needed comparator function */
  int (*compfunc)(const void *, const void *);

  /* pointer to the sort function, default is qsort */
  void (*sortfunc)(void *, size_t, size_t, int (*)(const void *, const void *))= qsort;

  /* process command-line parameters */
  /* the string here has the short option names, which are followed
     by a colon if they take a parameter, no colon if they are standalone */
  /* if an value exists after the option, it is specified in char *optarg */
  while ((ch = getopt_long(argc, argv, "f:d:rc:a:hv", 
         longopts, &option_index)) != -1) {
    switch (ch) {

    case 'f':
      /* the file name */
      strcpy(filename, optarg);
      break;

    case 'd':
      /* datatype */
      if (strcmp(optarg, "numeric") == 0) {
        datatype = NUMERIC;
      }
      else if (strcmp(optarg, "string") == 0) {
        datatype = STRING;
      }
      else if (strcmp(optarg, "ratio") == 0) {
        datatype = RATIO;
      }
      else {
        usage(argv[0]);
        return 1;
      }
      break;

    case 'c':
      /* criteria -- should do better error checking */
      if (strcmp(optarg, "length") == 0) {
        criterion = LENGTH;
      }
      else if (strcmp(optarg, "case-insensitive") == 0) {
        criterion = CASEINSENSITIVE;
      }
      else if (strcmp(optarg, "numerator") == 0) {
        criterion = NUMERATOR;
      }
      else {
        usage(argv[0]);
        return 1;
      }
      break;

    case 'a':
      /* algorithm -- should do better error checking */
      if (strcmp(optarg, "selection") == 0) {
        sortfunc = selsort;
      }
      else if (strcmp(optarg, "merge") == 0) {
        sortfunc = mergesort;
      }
      else if (strcmp(optarg, "quick") == 0) {
        sortfunc = qsort;
      }
      else if (strcmp(optarg, "heap") == 0) {
        sortfunc = heapsort;
      }
      else {
        usage(argv[0]);
        return 1;
      }
      break;

    case 'h':
      /* print help message */
      usage(argv[0]);
      return 0;

    case 'v':
      /* set verbose mode */
      verbose = 1;
      break;

    case 'r':
      /* set reverse mode */
      reverse = 1;
      break;

    default:
      usage(argv[0]);
      return 1;
    }
  }

  /* we have our command-line parameters parsed, let's get going */
  if (strcmp(filename, "") != 0) {
    infp = fopen(filename, "r");
    if (!infp) {
      fprintf(stderr, "%s: Could not open file %s for reading\n", argv[0], filename);
      return 1;
    }
  }
  else {
    infp = stdin;
  }

  /* read in actual data */
  int items = 0;
  switch (datatype) {

  case NUMERIC:
    while ((items < MAX_ELEMENTS) && fscanf(infp, "%d", &numbers[items]) == 1) {
      items++;
    }
    break;

  case STRING:
    while ((items < MAX_ELEMENTS) && fscanf(infp, "%s", word) == 1) {
      words[items] = (char *)malloc(strlen(word)+1);
      strcpy(words[items], word);
      while (fgetc(infp) != '\n');
      items++;
    }

    break;

  case RATIO:
    while ((items < MAX_ELEMENTS) && fscanf(infp, "%d/%d", &num, &den) == 2) {
      fractions[items] = create_ratio(num, den);
      items++;
    }

    break;
  }

  /* we have the data, let's sort! */
  switch (datatype) {

  case NUMERIC:
    if (reverse) {
      compfunc = intcmprev;
      if (criterion == LENGTH) {
        compfunc = intcmpdigitsrev;
      }
    }
    else {
      compfunc = intcmp;
      if (criterion == LENGTH) {
        compfunc = intcmpdigits;
      }
    }
    (*sortfunc)(numbers, items, sizeof(int), compfunc);
    break;

  case STRING:
    if (reverse) {
      compfunc = stringcmprev;
      if (criterion == LENGTH) {
        compfunc = stringcmplengthrev;
      }
      else if (criterion == CASEINSENSITIVE) {
        compfunc = stringcmpnocaserev;
      }
    }
    else {
      compfunc = stringcmp;
      if (criterion == LENGTH) {
        compfunc = stringcmplength;
      }
      else if (criterion == CASEINSENSITIVE) {
        compfunc = stringcmpnocase;
      }
    }

    
    (*sortfunc)(words, items, sizeof(char *), compfunc);
    break;

  case RATIO:
    if (reverse) {
      compfunc = ratiocmprev;
      if (criterion == NUMERATOR) {
        compfunc = ratiocmpnumrev;
      }
    }
    else {
      compfunc = ratiocmp;
      if (criterion == NUMERATOR) {
        compfunc = ratiocmpnum;
      }
    }

    (*sortfunc)(fractions, items, sizeof(ratio *), compfunc);
    break;
  }


  /* print the result */

  switch (datatype) {
  case NUMERIC:
    for (item = 0; item < items; item++) {
      printf("%d\n", numbers[item]);
    }
    break;

  case STRING:
    for (item = 0; item < items; item++) {
      printf("%s\n", words[item]);
    }
    break;

  case RATIO:
    for (item = 0; item < items; item++) {
      print_ratio(fractions[item]);
      printf("\n");
    }
    break;
  }

  // cleanup
  switch (datatype) {

  case STRING:
    for (item = 0; item < items; item++) {
      free(words[item]);
    }
    break;
  case RATIO:
    for (item = 0; item < items; item++) {
      free(fractions[item]);
    }
    break;
    case NUMERIC:
      // nothing to do
      break;
  }
  
  if (infp != stdin) fclose(infp);

  return 0;
}
