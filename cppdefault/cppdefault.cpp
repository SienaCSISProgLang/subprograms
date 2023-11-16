/*
  Example of optional parameters in a C++ program

  Jim Teresco
  The College of Saint Rose
  Programming Languages, CSC 433
  Fall 2012
*/

#include <iostream>

using namespace std;

// we provide default values for parameters should they be
// omitted in the call
int optional(int x, int y = 2, int z = 7) {

  cout << "x=" << x << ", y=" << y << ", z=" << z << endl;
  
}


int main(int argc, char *argv[]) {

  optional(1, 3, 8);
  optional(0, 9);
  optional(12);

  return 0;
}
