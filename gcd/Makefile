# Makefile for gcd program
#
# Jim Teresco, CS 433, The College of Saint Rose
#
# $Id$
#
CFILES=gcdmain.c gcd.c
OFILES=$(CFILES:.c=.o)
CC=gcc

gcdmain:	$(OFILES)
	$(CC) -o gcdmain $(OFILES)

clean::
	/bin/rm -f gcdmain $(OFILES)
