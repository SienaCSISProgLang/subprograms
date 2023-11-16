# Makefile for C++ example program
#
# Jim Teresco, CS 433, The College of Saint Rose
#
# $Id$
#
CXXFILES=callbyref.cpp
CFILES=cpointers.c
OXXFILES=$(CXXFILES:.cpp=.o)
OFILES=$(CFILES:.c=.o)
CXX=g++
CC=gcc

all:	callbyref cpointers

callbyref:	$(OXXFILES)
	$(CXX) -o callbyref $(OXXFILES)

cpointers:	$(OFILES)
	$(CC) -o cpointers $(OFILES)

clean::
	/bin/rm -f callbyref cpointers $(OFILES) $(OXXFILES)
