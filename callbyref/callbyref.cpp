/*
  Example of reference parameters in a C++ program

  Jim Teresco
  The College of Saint Rose
  Programming Languages, CSC 433
  Fall 2012
*/

#include <iostream>

using namespace std;

// x is passed by value, y by reference
void f(int x, int &y) {

  x = 23;
  y = 17;
}


int main(int argc, char *argv[]) {

  int a, b;
  a = 1;
  b = 2;

  cout << "Before: a=" << a << ", b=" << b << endl;

  f(a, b);

  cout << "After: a=" << a << ", b=" << b << endl;

  // this won't be allowed, since "2" can't be changed
  //f(1, 2);

  return 0;
}
