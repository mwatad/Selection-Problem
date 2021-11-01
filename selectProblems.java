
import java.util.concurrent.ThreadLocalRandom;

import javafx.util.Pair;

/**
*
* selectProblems
*
* @author Muhammad Watad
* @author Shadi Abu Shqara
* 
*/

public class selectProblems {
	
	/*
	 * Time complexity:
	 * Worst case: O(n^2)
	 * Best case: O(nlogn)
	 * Average case: O(nlogn)
	 */
	public Pair<Integer, Integer> selectRandQuickSort(int[] array, int k) {

		// Set the number of comparisons to 0
		Comparator.comparisons = 0;
		quickSort(array, 0, array.length - 1);

		return new Pair<Integer, Integer>(array[k - 1], (int) Comparator.comparisons);
		
	}
	
	/*
	 * Time complexity:
	 * Worst case: O(n^2)
	 * Best case: O(nlogn)
	 * Average case: O(nlogn)
	 * 
	 * Sorting sub-array of the input array in the range [l,r].
	 * The method uses hoare partition. 
	 */
	private static void quickSort(int[] array, int l, int r) {

		if (l < r) {
			//choosing pivot randomly 
			int i = randInt(l, r);
			// swap array[i] <-> array[r] if r is not the pivot's index
			if (i != r)
				swap(array, i, r);
			int p = hoarePartition(array, l, r);
			// Put the pivot in it's final place
			swap(array, p, r);
			// recursively sort right and left
			quickSort(array, l, p - 1);
			quickSort(array, p + 1, r);
		}

	}

	private static int hoarePartition(int[] array, int l, int r) {
		
		// choose the rightmost element to be the pivot
		int pivot = array[r];
		int i = l;
		int j = r - 1;

		while (true) {
			while (j >= l && Comparator.greaterThanOrEqual(array[j], pivot))
				j--;

			while (i <= r - 1 && Comparator.lessThanOrEqual(array[i], pivot))
				i++;

			// if i and j are stuck swap array[i] <-> array[j]
			if (i < j)
				swap(array, i, j);
			else
				return j + 1;
		}

	}

	private static void swap(int[] array, int i, int j) {

		int temp = array[i];
		array[i] = array[j];
		array[j] = temp;

	}

	// returns a random integer in the range [min, max]
	private static int randInt(int min, int max) {

		return ThreadLocalRandom.current().nextInt(min, max + 1);

	}

	/*
	 * Time complexity:
	 * Worst case: O(n^2)
	 * Best case: O(n)
	 * Average case: O(n^2)
	 */
	public Pair<Integer, Integer> selectInsertionSort(int[] array, int k) {

		// Set the number of comparisons to 0
		Comparator.comparisons = 0;
		insertionSort(array, 0, array.length - 1);

		return new Pair<Integer, Integer>(array[k - 1], (int) Comparator.comparisons);
	}

	/*
	 * Time complexity:
	 * Worst case: O(n^2)
	 * Best case: O(n)
	 * Average case: O(n^2)
	 */
	private static void insertionSort(int[] array, int l, int r) {

		for (int i = l + 1; i <= r; i++) {
			int temp = array[i];
			int j = i;
			
			while (j > l && Comparator.greaterThan(array[j - 1], temp)) {
				array[j] = array[j - 1];
				j--;
			}
			array[j] = temp;
		}

	}

	/*
	 * Time complexity:
	 * Worst case: O(n+klogn)
	 * Best case: O(n)
	 * Average case: O(n+klogn)
	 * 
	 * Select the k-th element by extracting the minimum
	 * k times from a minimum heap
	 */
	public Pair<Integer, Integer> selectHeap(int[] array, int k) {

		// Set the number of comparisons to 0
		Comparator.comparisons = 0;

		MinHeap minHeap = new MinHeap(array);

		int i = 1;
		while (i < k) {
			minHeap.extractMin();
			i++;
		}

		int value = minHeap.extractMin().getNum();

		return new Pair<Integer, Integer>(value, (int) Comparator.comparisons);

	}

	/*
	 * Time complexity:
	 * Worst case: O(n+klogk)
	 * Best case: O(n)
	 * Average case: O(n+klogk)
	 * 
	 * Select the k-th element of array using 2 minimum heaps.
	 * 
	 */
	public Pair<Integer, Integer> selectDoubleHeap(int[] array, int k) {

		// Set the number of comparisons to 0
		Comparator.comparisons = 0;

		MinHeap bigHeap = new MinHeap(array);
		MinHeap smallHeap = new MinHeap(array.length);

		// Insert the root of bigHeap to MinHeap
		// and save the root's index (pointer) 
		smallHeap.insert(bigHeap.getNode(0).getNum(), 0);
		int i = 1;
		while (i < k) {
			Node minimum = smallHeap.extractMin();
			i++;

			int leftIndex = bigHeap.left(minimum.getIndex());
			// If the extracted node has a left son
			// then add it to smallHeap
			if (leftIndex < bigHeap.size()) {
				int left = bigHeap.getNode(bigHeap.left(minimum.getIndex())).getNum();
				smallHeap.insert(left, leftIndex);
			}

			int rightIndex = bigHeap.right(minimum.getIndex());
			// If the extracted node has a right son
			// then add it to smallHeap
			if (rightIndex < bigHeap.size()) {
				int right = bigHeap.getNode(bigHeap.right(minimum.getIndex())).getNum();
				smallHeap.insert(right, rightIndex);
			}
		}

		int value = smallHeap.extractMin().getNum();

		return new Pair<Integer, Integer>(value, (int) Comparator.comparisons);

	}

	/*
	 * Time complexity:
	 * Worst case: O(n^2)
	 * Best case: O(n)
	 * Average case: O(n)
	 */
	public Pair<Integer, Integer> randQuickSelect(int[] array, int k) {
		// Set the number of comparisons to 0
		Comparator.comparisons = 0;
		int value = quickSelectAux(array, 0, array.length - 1, k - 1);

		return new Pair<Integer, Integer>(value, (int) Comparator.comparisons);
	}

	/*
	 * Time complexity:
	 * Worst case: O(n^2)
	 * Best case: O(n)
	 * Average case: O(n)
	 */
	private static int quickSelectAux(int[] array, int l, int r, int k) {

		// The array contains 1 element
		if (l == r)
			return array[l];

		// l < r
		
		// Choosing pivot randomly
		int i = randInt(l, r);
		swap(array, i, r);
		// Apply hoare partition
		int p = hoarePartition(array, l, r);
		// Put the pivot in it's final place
		swap(array, p, r);

		// The pivot is the k-th element
		if (p == k)
			return array[p];

		// The pivot is greater than k
		// continue recursively with the left sub-array
		if (p > k)
			return quickSelectAux(array, l, p - 1, k);
		
		// The pivot is smaller than k
		// continue recursively with the right sub-array
		else
			return quickSelectAux(array, p + 1, r, k);

	}

	
	/*
	 * Time complexity:
	 * Worst case: O(n)
	 * Best case: O(n)
	 * Average case: O(n)
	 */
	public Pair<Integer, Integer> medOfMedQuickSelect(int[] array, int k) {

		// Set the number of comparisons to 0
		Comparator.comparisons = 0;
		int value = medOfMedQuickSelectAux(array, 0, array.length - 1, k - 1);

		return new Pair<Integer, Integer>(value, (int) Comparator.comparisons);

	}

	/*
	 * Time complexity:
	 * Worst case: O(n)
	 * Best case: O(n)
	 * Average case: O(n)
	 */
	private static int medOfMedQuickSelectAux(int[] array, int l, int r, int k) {
		
		// number of elements in the range [l, r]
		int n = r - l + 1;
		
		// number of sub-arrays of size <= 5
		int[] medians = new int[(n + 4) / 5];

		// calculate medians for sub-arrays of size 5
		int i = 0;
		while (i < n / 5) {
			medians[i] = calculateMedian(array, l + i * 5, l + i * 5 + 4);
			i++;
		}

		// calculate the median for the the last sub-array
		// which has 1 to 4 elements only.
		if (i * 5 < n) {
			medians[i] = calculateMedian(array, l + i * 5, l + i * 5 + n % 5 - 1);
			i++;
		}

		int medOfMed;
		
		// There is only 1 median
		if (i == 1)
			medOfMed = medians[0];
		// calculate median of medians recursively
		else
			medOfMed = medOfMedQuickSelectAux(medians, 0, i - 1, i / 2);

		// Find the index of the median of medians in the range [l,r] 
		int medOfMedIndex = findIndex(array, medOfMed, l, r);
		swap(array, medOfMedIndex, r);
		int p = hoarePartition(array, l, r);
		swap(array, p, r);

		// the pivot is the k-th element
		if (p == k)
			return array[p];

		// the pivot is greater than the k-th element
		if (p > k)
			// continue recursivley on the left sub-array
			return medOfMedQuickSelectAux(array, l, p - 1, k);
		// the pivot is smaller than the k-th element
		else
			// continue recursivley on the right sub-array
			return medOfMedQuickSelectAux(array, p + 1, r, k);

	}

	private static int findIndex(int[] array, int val, int l, int r) {

		for (int i = l; i <= r; i++)
			if (Comparator.equal(array[i], val))
				return i;

		return -1;

	}

	private static int calculateMedian(int[] array, int l, int r) {

		insertionSort(array, l, r);
		return array[(l + r) / 2];

	}

}

class MinHeap {

	private Node[] elements;
	private int m;

	public MinHeap(int n) {
		elements = new Node[n];
		m = 0;
	}

	public MinHeap(int[] array) {
		elements = new Node[array.length];
		for (int e : array) {
			elements[m] = new Node(e, 0);
			m++;
		}

		// Heapify the internal nodes
		int k = array.length / 2 - 1;
		for (int i = k; i >= 0; i--)
			heapifyDown(i);
	}

	public void insert(int num, int index) {
		elements[m] = new Node(num, index);
		m++;
		heapifyUp(m - 1);
	}

	public Node extractMin() {

		Node min = elements[0];

		elements[0] = elements[m - 1];
		m--;
		heapifyDown(0);

		return min;
	}

	public Node getNode(int i) {
		return elements[i];
	}

	private void heapifyDown(int i) {

		int l = left(i);
		int r = right(i);
		int minimal = i;

		if (l < m && Comparator.lessThan(elements[l].getNum(), elements[minimal].getNum()))
			minimal = l;
		if (r < m && Comparator.lessThan(elements[r].getNum(), elements[minimal].getNum()))
			minimal = r;

		if (minimal != i) {
			swap(elements, minimal, i);
			heapifyDown(minimal);
		}
	}

	private void swap(Node[] elements, int i, int j) {

		Node temp = elements[i];
		elements[i] = elements[j];
		elements[j] = temp;

	}

	private void heapifyUp(int i) {

		int j = i;

		while (j > 0 && Comparator.lessThan(elements[j].getNum(), elements[parent(j)].getNum())) {
			swap(elements, j, parent(j));
			j = parent(j);
		}

	}

	private int parent(int i) {

		if (i % 2 == 1)
			return i / 2;

		return (i - 1) / 2;

	}

	public int right(int i) {
		return 2 * i + 2;
	}

	public int left(int i) {
		return 2 * i + 1;
	}

	public int size() {
		return m;
	}

}

class Node {

	private int num;
	private int index; //index of a node

	public Node(int num, int index) {
		this.num = num;
		this.index = index;
	}

	public int getNum() {
		return num;
	}

	public void setNum(int num) {
		this.num = num;
	}

	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}

}

/*
 * This class is used to count the number of comparisons
 * We replaced the operations <  , <= , > , >= , == , with
 * operations that updates a counter.
 */
class Comparator {

	public static long comparisons = 0;

	public static boolean lessThan(int x, int y) {
		comparisons++;
		return x < y;
	}

	public static boolean lessThanOrEqual(int x, int y) {
		comparisons++;
		return x <= y;
	}

	public static boolean greaterThan(int x, int y) {
		comparisons++;
		return x > y;
	}

	public static boolean greaterThanOrEqual(int x, int y) {
		comparisons++;
		return x >= y;
	}

	public static boolean equal(int x, int y) {
		comparisons++;
		return x == y;
	}

}
