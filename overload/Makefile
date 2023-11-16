# Makefile for C++ example program
#
# Jim Teresco, CS 433, The College of Saint Rose
#
# $Id$
#
CXXFILES=overload.cpp
OFILES=$(CXXFILES:.cpp=.o)
CXX=g++

overload:	$(OFILES)
	$(CXX) -o overload $(OFILES)

clean::
	/bin/rm -f overload $(OFILES)
