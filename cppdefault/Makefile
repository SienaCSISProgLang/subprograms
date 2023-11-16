# Makefile for C++ example program
#
# Jim Teresco, CS 433, The College of Saint Rose
#
# $Id$
#
CXXFILES=cppdefault.cpp
OFILES=$(CXXFILES:.cpp=.o)
CXX=g++

cppdefault:	$(OFILES)
	$(CXX) -o cppdefault $(OFILES)

clean::
	/bin/rm -f cppdefault $(OFILES)
