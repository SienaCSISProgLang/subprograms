/*
  SimpleLinkedListObject.java
  
  Implementation of a very simple singly-linked list structure using 
  object type (not generics)
  
  Jim Teresco, The College of Saint Rose, CSC 433, Fall 2014

  Based on previous implementations from
  CSC 136, Williams College, COMSC 211, Mount Holyoke College, and CSC
  501 Fall 2013 and CSC 523, Summer 2014, The College of Saint Rose
  
  $Id$
*/

import java.util.Iterator;

/*
The list node: a class which is not public -- only the other classes
within this file will be able to construct one.
 */

class SimpleListNodeObject {
    /** Object to be stored in this node */
    protected Object value;
    /** reference to next node, null if this is the last one */
    protected SimpleListNodeObject next;

    /**
    construct a new list object

    @param value value to be stored in this list node
    @param next reference to next node in the list
     */
    public SimpleListNodeObject(Object value, SimpleListNodeObject next) {
        this.value = value;
        this.next = next;
    }

    /**
    get data value from this list object

    @return data object stored
     */
    public Object value() {

        return value;
    }

    /**
    get next reference from this list object

    @return reference to next list object
    @return null if this is the last object in the list
     */
    public SimpleListNodeObject next() {

        return next;
    }

    /**
    set the data value of this object to a new value.  Old value is returned

    @param value new data value to be stored
    @return old stored data value
     */
    public Object setValue(Object value) {
        Object prev = this.value;
        this.value = value;
        return prev;
    }

    /**
    set the next reference of this object to a new value.  
    Old next reference is returned

    @param next new next reference to be stored
    @return old stored next reference
     */
    public SimpleListNodeObject setNext(SimpleListNodeObject next) {
        SimpleListNodeObject prev = this.next;
        this.next = next;
        return prev;
    }
}

/**
Simple singly-linked list structure's iterator.

 */
class SimpleListIteratorObject implements java.util.Iterator {
    /** Reference to currently considered element within list */
    protected SimpleListNodeObject current;
    /** The head of list for reset */
    protected SimpleListNodeObject head;

    /**
    Construct an iterator that traverses list beginning at t.

    @post returns an iterator that traverses a linked list

    @param t The first element of list to be traversed.
     */
    public SimpleListIteratorObject(SimpleListNodeObject t) {
        head = t;
        current = head;
    }

    /**
    Determine if the iteration is finished.

    @post returns true if there is more structure to be viewed:
    i.e., if value (next) can return a useful value.

    @return True if the iterator has more elements to be considered.
     */
    public boolean hasNext() {
        return current != null;
    }

    /**
    Return current value and increment Iterator.

    @pre traversal has more elements
    @post returns current value and increments iterator

    @return The current value, before increment.
     */
    public Object next() {
        Object temp = current.value();
        current = current.next();
        return temp;
    }

    /**
    Remove operation not supported by this Iterator
     */
    public void remove() {

        throw new UnsupportedOperationException();
    }
}

/**
Generic singly-linked list structure.  This has only simple nodes,
no tail pointer is included.  This uses Java 1.5 generics.

 */
public class SimpleLinkedListObject implements Iterable {

    /** reference to the node that contains the first value in the list */
    protected SimpleListNodeObject head;

    /**
    construct a new, empty list
     */
    public SimpleLinkedListObject() {
        clear();
    }

    /**
    add a new item as the pos'th object in the list

    @param pos desired position of the object in the list
    @param obj data object to insert
    @pre the list must have between 0 and pos-1 entries
    @post the list is modified to have obj at position pos, subseqent items slide down
     */
    public void add(int pos, Object obj) {

	// negative positions not allowed
	if (pos < 0) {
	    throw new IndexOutOfBoundsException("Attempt to add at negative position " + pos);
        }

        // if we are adding to an empty list it must be at position 0
        if ((head == null) && (pos != 0)) {
	    throw new IndexOutOfBoundsException("Attempt to add at position " + pos + " in empty list");
        }

        // are we adding at the front?  note that if head is null, this
	// inserts a first element, if not, it inserts at the head of
	// an existing list
        if (pos == 0) {
            head = new SimpleListNodeObject(obj, head);
            return;
        }

        // we are adding somewhere else, find entry after which we will 
        // insert our item
        int i = 0;
        SimpleListNodeObject finger = head;
        while (i < pos-1) {
            i++;
            finger = finger.next();
            if (finger == null) {
                throw new IndexOutOfBoundsException("Attempt to add at position " + pos + " in list of size " + i);
            }
        }
        // finger points at the node after which we want to add
        // so create the new object with finger's next as its next
        // and set finger's next to the new node.
        // note that this also works for the case when we are adding
        // to the end
        finger.setNext(new SimpleListNodeObject(obj, finger.next()));
    }

    /**
       add at 0 is the default kind of add, and it simply calls the
       general purpose add

       @param the item to be added to the start of the list
    */
    public void add(Object obj) {

	add(0, obj);
    }

    /**
    retrieve an item in the list

    @param pos position at which to find the element
    @pre the list has at least pos+1 elements
    @return the item at position pos
     */
    public Object get(int pos) {

	// negative positions not allowed
	if (pos < 0) {
	    throw new IndexOutOfBoundsException("Attempt to get from a negative position " + pos);
        }

        SimpleListNodeObject finger = head;
        int i = 0;

        if (head == null) {
            throw new IndexOutOfBoundsException("Attempt to get from an empty list");
        }

        while (i < pos) {
            i++;
            finger = finger.next();
            if (finger == null) {
                throw new IndexOutOfBoundsException("Attempt to get element " + pos + " from a " + i + " element list");
            }
        }
        return finger.value();
    }

    /**
    set an item in the list by position

    @param pos position at which to replace an element
    @pre the list has at least pos+1 elements
    @return the previous item at position pos
     */
    public Object set(int pos, Object obj) {

	// negative positions not allowed
	if (pos < 0) {
	    throw new IndexOutOfBoundsException("Attempt to set at negative position " + pos);
        }

        SimpleListNodeObject finger = head;
        int i = 0;

        if (head == null) {
            throw new IndexOutOfBoundsException("Attempt to set in an empty list");
        }

        while (i < pos) {
            i++;
            finger = finger.next();
            if (finger == null) {
                throw new IndexOutOfBoundsException("Attempt to set element " + pos + " from a " + i + " element list");
            }
        }
        return finger.setValue(obj);    
    }

    /**
    search for an item in the list, return whether it is found

    @param obj object to search for in list
    @return true if the item is found
    @return false if the items is not in the list
     */
    public boolean contains(Object obj) {

        // easy when the list is empty
        if (head == null) return false;

        // otherwise look for it
        SimpleListNodeObject finger = head;
        while (finger != null) {
            if (finger.value().equals(obj)) return true;
            finger = finger.next();
        }
        return false;
    }

    /**
    return the number of elements in the list

    @return the number of elements in the list
     */
    public int size() {
        int count = 0;

        SimpleListNodeObject finger = head;
        while (finger != null) {
            count++;
            finger = finger.next();
        }

        return count;
    }

    /**
    remove an item from the list, by position

    @param pos the position whose element is to be removed
    @return the object removed
     */
    public Object remove(int pos) {

	// negative positions not allowed
	if (pos < 0) {
	    throw new IndexOutOfBoundsException("Attempt to remove from negative position " + pos);
        }

        if (head == null) {
            throw new IndexOutOfBoundsException("Attempt to remove from an empty list");
        }

        // check for removal of the first item in the list
        // this works for all cases including the case of a 
	// one-element list, as head gets set to null
        if (pos == 0) {
            Object retval = head.value();
            head = head.next();
            return retval;
        }

        // remove an item at a non-first position
        SimpleListNodeObject finger = head;
        int count = 0;
        // find the item before the one we want to remove
        while (count < pos-1) {
            count++;
            finger = finger.next();
            if (finger == null) {
                throw new IndexOutOfBoundsException("Attempt to remove element at index " + pos + " in a list with " + count + " elements");
            }
        }
        // finger is pointing to item pos-1
        // make sure there is something at pos
        if (finger.next() == null) {
            throw new IndexOutOfBoundsException("Attempt to remove element at index " + pos + " in a list with " + count + " elements");
        }
        Object retval = finger.next().value();
        finger.setNext(finger.next().next());
        return retval;
    }

    /**
    remove all items from the list

    @post all items are removed from the list
     */
    public void clear() {

        head = null;
    }

    /**
    create iterator for the list

    @return iterator over SimpleLinkedList
     */
    public Iterator iterator() {

        return new SimpleListIteratorObject(head);
    }

    /**
    generate a string representation of a linked list 

    @return a representation of the list items as a String
     */
    public String toString() {

        StringBuffer s = new StringBuffer();
        s.append("<SimpleLinkedListObject:");
        Iterator li = iterator();
        while (li.hasNext())
        {
            s.append(" "+li.next());
        }
        s.append(">");
        return s.toString();
    }

    public static void main(String[] args) {

	// build a list and do some things with it
        SimpleLinkedListObject list = new SimpleLinkedListObject();

        list.add(0, 1);
        list.add(1, 2);
        list.add(2, 3);
        list.add(4);
        list.add(2, 5);
        list.set(3, 6);

        System.out.println(list);

        System.out.println("Size: " + list.size());
        System.out.println("Item at position 2: " + list.get(2));
        System.out.println("contains(3): " + list.contains(3));
        System.out.println("contains(4): " + list.contains(4));
        System.out.println("contains(5): " + list.contains(5));
        System.out.println("contains(23): " + list.contains(23));

        System.out.println("remove(3) removed : " + list.remove(3));
        System.out.println(list);

        System.out.println("remove(0) removed : " + list.remove(0));
        System.out.println(list);

        int first = (Integer)list.get(0);
        int second = (Integer)list.get(1);

        // look, a for each loop!  we can do this because our class
	// implements the Iterable interface.
        int sum = 0;
        for (Object i : list) {
            sum += (Integer)i;
        }
        System.out.println("Sum of entries: " + sum);
    }
}
