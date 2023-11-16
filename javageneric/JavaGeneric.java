/*
  Example of generic methods in Java

  Jim Teresco
  The College of Saint Rose
  Programming Languages, CSC 433
  Fall 2012

  Based on examples in Sebesta.
*/

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;

public class JavaGeneric {
    
    // a method that takes a type parameter and a parameter of that type
    public static <T> void printInfo(T obj) {

	// not much we can do with an object of an arbitrary type,
	// but Java's Object class does provide a few
	// this calls toString
	System.out.println("printInfo toString: " + obj);
	// hash code
	System.out.println("printInfo hashCode: " + obj.hashCode());
	// print the object's class
	System.out.println("printInfo getClass: " + obj.getClass());

    }

    // this one places the restriction that the type parameter
    // must be one that extends the Comparable class, meaning
    // it must provide a compareTo method
    public static <T extends Comparable> T maxOfArray(T[] a) {

	// handle empty array case
	if (a.length == 0) return null;

	// we can declare locals of the type parameter
	T winner = a[0];

	for (int i=1; i<a.length; i++) {
	    if (winner.compareTo(a[i]) < 0) {
		winner = a[i];
	    }
	}

	return winner;
    }

    // this one only works on ArrayLists of some type -- a wildcard type
    // is used to specify this
    public static void traverseArrayList(ArrayList<?> al) {

	for (Object item : al) {
	    System.out.println("traverseArrayList: item is " + item);
	}
    }

    // this one is more generic - works on any class that implements the
    // Collection interface
    public static void traverseCollection(Collection<?> c) {

	for (Object item : c) {
	    System.out.println("traverseCollection: item is " + item);
	}
    }


    public static void main(String args[]) {

	// try out printInfo
	printInfo(new Integer(17));
	printInfo("Some String!");
	printInfo(new ArrayList());

	// try out maxOfArray
	Integer[] a1 = {1, 3, 9, 12, 9, 5, 0};
	System.out.println("maxOfArray(a1) = " + maxOfArray(a1));
	String[] a2 = new String[6];
	a2[0] = "Maine";
	a2[1] = "New Hampshire";
	a2[2] = "Vermont";
	a2[3] = "Massachusetts";
	a2[4] = "Rhode Island";
	a2[5] = "Connecticut";
	System.out.println("maxOfArray(a2) = " + maxOfArray(a2));

	// we wouldn't be able to use it on some other arrays
	// such as an array of LinkedLists (if we had one of
	// those for some reason
	LinkedList[] a3 = new LinkedList[5];
	a3[0] = new LinkedList();
	a3[1] = new LinkedList();
	a3[2] = new LinkedList();
	a3[3] = new LinkedList();
	a3[4] = new LinkedList();
	// would cause error:
	//System.out.println("maxOfArray(a3) = " + maxOfArray(a3));

	// try out the traverseArrayList method
	ArrayList<Integer> a4 = new ArrayList<Integer>();
	a4.add(1);
	a4.add(800);
	a4.add(123);
	a4.add(5678);
	traverseArrayList(a4);

	ArrayList<String> a5 = new ArrayList<String>();
	a5.add("Programming");
	a5.add("Languages");
	a5.add("Rules");
	traverseArrayList(a5);

	// now the second as a Collection instead (since ArrayList implements
	// Collection
	traverseCollection(a5);

	LinkedList<Double> a6 = new LinkedList<Double>();
	a6.add(0.0);
	a6.add(Math.PI);
	a6.add(Math.sqrt(2.0));
	traverseCollection(a6);


    }	
}
