# Makefile for varargs example
#
# Jim Teresco, CS 433, The College of Saint Rose
#
# $Id$
#
CFILES=varargs.c
OFILES=$(CFILES:.c=.o)
CC=gcc -std=c99

varargs:	varargs.o
	$(CC) -o varargs varargs.o

clean::
	/bin/rm -f varargs $(OFILES)
