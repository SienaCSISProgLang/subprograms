# Makefile for jensens example program
#
# Jim Teresco, CS 433, The College of Saint Rose
#
# $Id$
#
CFILES=jensens.c
OFILES=$(CFILES:.c=.o)
CC=gcc -std=c99

jensens:	$(OFILES)
	$(CC) -o jensens $(OFILES)

clean::
	/bin/rm -f jensens $(OFILES)
