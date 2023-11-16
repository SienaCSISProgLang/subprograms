/*
  Example of a generic function in a C++ program

  Jim Teresco
  The College of Saint Rose
  Programming Languages, CSC 433
  Fall 2012

  "max" function based on Sebesta.
*/

#include <iostream>
#include <string>

using namespace std;

// a generic "max" function -- should work on any datatype which has
// instances that can be compared with >
template <class Type>
  Type mymax(Type first, Type second) {
     return first > second ? first : second;
  }


// and we try a few out
int main(int argc, char *argv[]) {

  cout << "mymax(1, 3) = " << mymax(1, 3) << endl;

  cout << "mymax(5.6, 3.14) = " << mymax(5.6, 3.14) << endl;

  cout << "mymax('A', '3') = " << mymax('A', '3') << endl;

  // these compare the pointers, not the strings...
  cout << "mymax(\"hello\", \"world\") = " << mymax("hello", "world") << endl;
  cout << "mymax(\"goodbye\", \"apple\") = " << mymax("goodbye", "apple") << endl;
  string first = "alice", second = "bob", third = "carol", fourth = "dilbert";

  // these should compare the strings, using operator< from the string class
  cout << "mymax(" << first << ", " << second << ") = " << mymax(first,second) << endl;
  cout << "mymax(" << fourth << ", " << third << ") = " << mymax(fourth,third) << endl;

  return 0;
}
