/*
  Variable argument list length example in Java

  Jim Teresco
  The College of Saint Rose
  Programming Languages, CSC 433
  Fall 2014
*/

public class Varargs {

    // variable length parameter list of ints to add up -- no need for a
    // sentinel here, since the variable list comes to the method as an 
    // array.
    // the varargs parameter must be last in the parameter list.
    public static int sumNumbers(int... numbers) {

	// now loop over the arguments and accumulate the sum
	int sum = 0;
	for (int number : numbers)  {
	    sum += number;
	}
	
	return sum;
    }
    
    public static void main(String args[]) {

	System.out.println("First: " + sumNumbers(1, 2, 3, -1));
	System.out.println("Second: " + sumNumbers(82, 12, 21, 211, 98, 7, 12, 87, -23));
	System.out.println("Fourth: " + sumNumbers(1, 2, 3, 4, 5, 6, 7, 8, 9, 10));
	// what if we give it the wrong parameters?
	//System.out.println("Fifth: " + sumNumbers(1, 8.2, "hi there!", 87L, -1));
    }
}

