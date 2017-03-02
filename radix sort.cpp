#include <iostream>

using namespace std;

int const bitsSortedOn = 3;

// optimized standard radix sort without copy back
int radixSortWithoutCopyBack(int *arr, int size, int bitsSortedOn) {
	int const buckets = 1 << bitsSortedOn;
	int const k = buckets-1;
	int freq[buckets];
	for (int i=0; i<buckets; i++) freq[i] = 0;
	int *tmp = new int[size];
	int shift = 0;

	while (shift < 32) {
		for (int i=0; i<size; i++) { // count frequencies
			freq[(arr[i] >> shift) & k]++; }
		for (int i=1; i<buckets; i++) { // sum frequencies 
			freq[i] += freq[i-1]; }
		for (int i=size-1; i >= 0; i--) { // move nodes correct loc in tmp array
			tmp[--freq[(arr[i] >> shift) & k]] = arr[i]; }
		for (int i=0; i<buckets; i++) { // set frequencies to 0
			freq[i] = 0; }

		shift += bitsSortedOn;
		if (!(shift < 32)) { break; }

		for (int i=0; i<size; i++) { // count frequencies
			freq[(tmp[i] >> shift) & k]++; }
		for (int i=1; i<buckets; i++) { // sum frequencies
			freq[i] += freq[i-1]; }
		for (int i=size-1; i >= 0; i--) { // move nodes to correct loc in arr
			arr[--freq[(tmp[i] >> shift) & k]] = tmp[i]; }
		for (int i=0; i<buckets; i++) { // set frequencies to 0
			freq[i] = 0; }
		shift += bitsSortedOn;

	}
	if (32 % bitsSortedOn != 0) {
		for (int i=0; i<size; i++) {
			arr[i] = tmp[i];
		}
	}
	delete[] tmp;

	return 0;
}


int radixSort(int *arr, int size, int bitsSortedOn) {
	int const buckets = 1 << bitsSortedOn;
	int const k = buckets-1;
	int freq[buckets]; 
	for (int i=0; i<buckets; i++) freq[i] = 0;
	int *tmp = new int[size];
	int shift = 0;

	while (shift < 32) {
		for (int i=0; i<size; i++) { // count frequencies
			freq[(arr[i] >> shift) & k]++; }
		for (int i=1; i<buckets; i++) { // sum frequencies 
			freq[i] += freq[i-1]; }
		for (int i=size-1; i >= 0; i--) { // move nodes correct loc in tmp array
			tmp[--freq[(arr[i] >> shift) & k]] = arr[i]; }
		for (int i=0; i<buckets; i++) { // set frequencies to 0
			freq[i] = 0; }
		for (int i=0; i<size; i++) { // copy from tmp back to arr
			arr[i] = tmp[i]; }
		shift += bitsSortedOn;
	}

	delete[] tmp;
	
	return 0;
}

/*
 * Idea from algorithm 1 in article
 * Allocate buckets for each possible value for the sorted in bits
 * Uses much more memory when sorting on more bits
 */
int radixSortWoCountingFreq(int *arr, int size, int bitsSortedOn) {
	int const buckets = 1 << bitsSortedOn;
	int const k = buckets-1;
	int freq[buckets];
	for (int i=0; i<buckets; i++) freq[i] = 0;
	// can't use new since buckets and size aren't known and compile time
	int *tmp = (int *) malloc(sizeof(int) * size * buckets);
	int shift = 0;

	while (shift < 32) {
		for (int i=0; i < size; i++) {
			int bucket = (arr[i] >> shift) & k;
			tmp[freq[bucket] + size * bucket] = arr[i];
			freq[bucket]++;
		}
		// Go through each bucket and copy all its content back to arr
		int elem = size-1;
		for (int bucket=buckets-1; bucket>=0; bucket--) {
			while (freq[bucket] != 0) {
				arr[elem--] = tmp[bucket*size + --freq[bucket]];
			}
		}
		for (int i=0; i<buckets; i++) { // set frequencies to 0
			freq[i] = 0; }
		shift += bitsSortedOn;
	}

	free(tmp);

	return 0;
}

bool isSorted(int *arr, int size) {
	for (int i=0, end = size-1; i<end; i++) {
		if (arr[i] > arr[i+1]) {
			return false;
		}
	}
	return true;
}

void testBitsSortedOn(int N, int reps) {
	int *array = new int[N];
	for (int i = 0; i<N; i++) {
		array[i] = rand();
	}
	int times[17] = {0}; // larger than 16 doesn't make really make sense
	for (int k=0; k<=reps;k++) {
		for (int i=1; i<=16;i++) {
			for (int j = 0; j<N; j++) {
				array[j] = rand();
			}
			clock_t t = clock();
			radixSort(array, N, i);
			times[i] += clock() - t;
		}
	}
	for (int i = 1; i<17; i++) {
		cout << (times[i]/reps) << " when sorting on " << i << " bits" << endl;
	}
	delete[] array;
}

int main() { }
