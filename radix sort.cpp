#include <iostream>

using namespace std;

int const bitsSortedOn = 3;

// optimized standard radix sort without copy back
int radixSortWithoutCopyBack(int *arr, int size, int bitsSortedOn) {
	int const k = (1 << bitsSortedOn)-1;
	int const buckets = k+1;
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
	int const k = (1 << bitsSortedOn)-1;
	int const buckets = k+1;
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
	int times[17] = {0};
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

int main() {}
