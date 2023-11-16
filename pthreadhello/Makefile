# Makefile for pthread hello world program
#
# Jim Teresco, CS 338, Williams College
# CS 341, Mount Holyoke College
# CS 400, Siena College
#
# $Id$
#
CFILES=pthreadhello.c
OFILES=$(CFILES:.c=.o)
# On Solaris, cc -mt will build an executable with threading support
CC=cc -mt

pthreadhello:	$(OFILES)
	$(CC) -o pthreadhello $(OFILES) -lpthread

clean::
	/bin/rm -f pthreadhello $(OFILES)
