/*
  Example of overloaded functions in a C++ program

  Jim Teresco
  The College of Saint Rose
  Programming Languages, CSC 433
  Fall 2012
*/

#include <iostream>

using namespace std;

// we have 3 functions f, and as long as they all have different
// protocols, C++ can choose the right one based on the call
void f(int x, int y) {

  cout << "It's the first f" << endl;
}


void f(int x) {

  cout << "It's the second f" << endl;
}


void f(char *str) {

  cout << "It's the third f" << endl;
}


int main(int argc, char *argv[]) {

  f(1, 2);
  f(3);
  f("hi there!");

  return 0;
}
