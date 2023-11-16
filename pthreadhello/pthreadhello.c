/*
  Simple pthread program

  Jim Teresco

  Fall 2000, Modified Fall 2002, Spring 2010, Spring 2012

*/

/* include some stuff we need */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
/* including the pthread header file */
#include <pthread.h>

/* This is going to be our thread function */
void *thread(void *args) {

  int *val=(int *)args;
  /* we just identify ourselves by the value sent in when the thread
     was created (which happens below in the main program */
  printf("This is thread with val %d\n", *val);

  /* when a thread is going to exit, it's a good idea to call pthread_exit() */
  pthread_exit(0);
}

/* The main program starts here */
int main(int argc, char *argv[]) {

  /* we will create two threads below.  We declare two pthread_t
     variables to remember their threadids. */
  pthread_t child1, child2;
  /* allocate some space to hold the values we'll pass to the thread
     functions */
  int vals[] = {1, 2};
  int rc;


  /* The main thread announces its existence */
  printf("Hello from the main thread\n");

  /* We create the child threads.  The third argument specifies the
     function to call when the new thread is created */
  rc=pthread_create(&child1, NULL, thread, (void *)&(vals[0]));
  /* Check to make sure it worked */
  if (rc != 0) {
    fprintf(stderr,"Could not create child thread 1\n");
  }

  /* And do it again */
  rc=pthread_create(&child2, NULL, thread, (void *)&(vals[1]));
  if (rc != 0) {
    fprintf(stderr,"Could not create child thread 2\n");
  }

  /* Wait for each child thread to call pthread_exit() */
  pthread_join(child1,NULL);
  pthread_join(child2,NULL);
 
  /* Print one more message.. */
  printf("Goodbye from the main thread\n");

  /* and be done. */
  return 0;
}
