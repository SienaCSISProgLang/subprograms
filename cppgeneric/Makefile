# Makefile for C++ example program
#
# Jim Teresco, CS 433, The College of Saint Rose
#
# $Id$
#
CXXFILES=cppgeneric.cpp
OFILES=$(CXXFILES:.cpp=.o)
CXX=g++

cppgeneric:	$(OFILES)
	$(CXX) -o cppgeneric $(OFILES)

clean::
	/bin/rm -f cppgeneric $(OFILES)
