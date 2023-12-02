# Makefile for C sorting demo program
#
# Jim Teresco
#
# CSC 381, The College of Saint Rose
# Updated for CSIS-340, Siena College
#
CFILES=ratio.c gcd.c csorts.c
OFILES=$(CFILES:.c=.o)
CC=gcc -g -Wall

csorts:	$(OFILES)
	$(CC) -o csorts $(OFILES)

clean::
	/bin/rm -f csorts $(OFILES)
