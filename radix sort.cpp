#include <iostream>
#include <bitset>
//#include <emmintrin.h>
#include <xmmintrin.h>

using namespace std;

int const bitsSortedOn = 3;

void printBinaryArray(int *arr, int n) {
	for (int i=0; i<n; i++) {
		cout << bitset<32>(arr[i]) << endl;
	}
}

void printBucketArray(int *matrix, int n, int *freq, int bitsSortedOn) {
	int buckets = 1 << bitsSortedOn;
	for (int i = 0; i < buckets; i++) {
		int bucketOffset = i * n, j = 0;
		cout << "printing contents of bucket " << i << endl;
		while (j < freq[i]) {
			cout << "    " << hex << matrix[bucketOffset + j++] << endl;
		}
	}
	cout << endl;
}
void printHexArray(int *arr, int n) {
	for (int i=0; i<n; i++) {
		cout << hex << arr[i] << endl;
	}
}

void printArray(int *arr, int n) {
	for (int i=0; i<n; i++) {
		cout << arr[i] << endl;
	}
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
			int index = --freq[(arr[i] >> shift) & k];
			// tmp[index] = arr[i]; 
			_mm_stream_si32(&tmp[index], arr[i]);
		}
		for (int i=0; i<buckets; i++) { // set frequencies to 0
			freq[i] = 0; }
		for (int i=0; i<size; i++) { // copy from tmp back to arr
			arr[i] = tmp[i]; }
		shift += bitsSortedOn;
	}

	delete[] tmp;

	return 0;
}

//
int radixSortFreqMatrix(int *arr, int size, int bitsSortedOn) {
	int const buckets = 1 << bitsSortedOn;
	int const k = buckets-1;
	int const passes = 32 / bitsSortedOn + 1;
	// better to store on heap with optimized zero assignment?
	int freq[passes][buckets];
	for (int i=0; i<passes; i++)
		for (int j=0; j<buckets; j++)
			freq[i][j] = 0;
	int *tmp = new int[size];
	int shift = 0, iteration = 0;

	while (shift < 32) {
		for (int i=0; i<size; i++) { // count frequencies
			freq[iteration][(arr[i] >> shift) & k]++; }
		for (int i=1; i<buckets; i++) { // sum frequencies
			freq[iteration][i] += freq[iteration][i-1]; }
		for (int i=size-1; i >= 0; i--) { // move nodes correct loc in tmp array
			tmp[--freq[iteration][(arr[i] >> shift) & k]] = arr[i]; }
		for (int i=0; i<size; i++) { // copy from tmp back to arr
			arr[i] = tmp[i]; }
		shift += bitsSortedOn;
		iteration++;
	}

	delete[] tmp;

	return 0;
}

// optimized standard radix sort without copy back
int radixSortWoCopyBack(int *arr, int size, int bitsSortedOn) {
	int const buckets = 1 << bitsSortedOn;
	int const k = buckets-1;
	int freq[buckets];
	for (int i=0; i<buckets; i++) freq[i] = 0;
	int *tmp = new int[size];
	int shift = 0;
	bool exitedEarly = false;

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
		if (!(shift < 32)) { exitedEarly = true; break; }

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
	if (exitedEarly) {
		for (int i=0; i<size; i++) {
			arr[i] = tmp[i];
		}
	}
	delete[] tmp;

	return 0;
}

int radixSortWoCopyBackFreqMatrix(int *arr, int size, int bitsSortedOn) {
	int const buckets = 1 << bitsSortedOn;
	int const k = buckets-1;
	int const passes = 32 / bitsSortedOn + 1;
	int freq[passes][buckets];
	for (int i=0; i<passes; i++)
		for (int j=0; j<buckets; j++)
			freq[i][j] = 0;
	int *tmp = new int[size];
	int shift = 0, iteration = 0;
	bool exitedEarly = false;

	while (shift < 32) {
		for (int i=0; i<size; i++) { // count frequencies
			freq[iteration][(arr[i] >> shift) & k]++; }
		for (int i=1; i<buckets; i++) { // sum frequencies
			freq[iteration][i] += freq[iteration][i-1]; }
		for (int i=size-1; i >= 0; i--) { // move nodes correct loc in tmp array
			tmp[--freq[iteration][(arr[i] >> shift) & k]] = arr[i]; }

		iteration++;
		shift += bitsSortedOn;
		if (!(shift < 32)) { exitedEarly = true; break; }

		for (int i=0; i<size; i++) { // count frequencies
			freq[iteration][(tmp[i] >> shift) & k]++; }
		for (int i=1; i<buckets; i++) { // sum frequencies
			freq[iteration][i] += freq[iteration][i-1]; }
		for (int i=size-1; i >= 0; i--) { // move nodes to correct loc in arr
			arr[--freq[iteration][(tmp[i] >> shift) & k]] = tmp[i]; }
		shift += bitsSortedOn;
		iteration++;
	}
	if (exitedEarly) {
		for (int i=0; i<size; i++) {
			arr[i] = tmp[i];
		}
	}
	delete[] tmp;

	return 0;
}

/*
 * Idea from algorithm 1 in article
 * Allocate buckets for each possible value for the sorted in bits
 * Uses much more memory when sorting on more bits
 *
 * segfaults on my machine, macbook pro 2011 w 16gb ram, when sorting this many ints
 * 15 bits on 100.000 numbers
 * 12 bits on 1.000.000 numbers
 * 8 bits on 10.000.000 numbers
 * 5 bits on 100.000.000 numbers
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
			long bucket = (arr[i] >> shift) & k;
			tmp[freq[bucket] + size * bucket] = arr[i];
			freq[bucket]++;
		}
		// Go through each bucket and copy all its content back to arr
		int elem = size-1;
		for (long bucket=buckets-1; bucket>=0; bucket--) {
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

// Same behaviour as the version above
int radixSortWoCountingFreqWFreqMatrix(int *arr, int size, int bitsSortedOn) {
	int const buckets = 1 << bitsSortedOn;
	int const k = buckets-1;
	int const passes = 32 / bitsSortedOn + 1;
	int freq[passes][buckets];
	for (int i=0; i<passes; i++)
		for (int j=0; j<buckets; j++)
			freq[i][j] = 0;
	// can't use new since buckets and size aren't known and compile time
	int *tmp = (int *) malloc(sizeof(int) * size * buckets);
	int shift = 0, iteration = 0;

	while (shift < 32) {
		for (int i=0; i < size; i++) {
			long bucket = (arr[i] >> shift) & k;
			tmp[freq[iteration][bucket] + size * bucket] = arr[i];
			freq[iteration][bucket]++;
		}
		// Go through each bucket and copy all its content back to arr
		int elem = size-1;
		for (long bucket=buckets-1; bucket>=0; bucket--) {
			while (freq[iteration][bucket] != 0) {
				arr[elem--] = tmp[bucket*size + --freq[iteration][bucket]];
			}
		}
		shift += bitsSortedOn;
		iteration++;
	}

	free(tmp);

	return 0;
}

/*
 * With two matrixes to sort back and forth between
 * freq matrix can be change to two vectors
 * can't sort up to 1.000.000 without problems
 * above that it seems to have same behaviour as with only 1 tmp matrix
*/
int radixSortWoCountingFreqW2Tmps(int *arr, int size, int bitsSortedOn) {
	int const buckets = 1 << bitsSortedOn;
	int const k = buckets-1;
	int const passes = 32 / bitsSortedOn + 1;
	int freq[passes][buckets];
	for (int i=0; i<passes; i++)
		for (int j=0; j<buckets; j++)
			freq[i][j] = 0;
	// can't use new since buckets and size aren't known and compile time
	int *tmp1 = (int *) malloc(sizeof(int) * size * buckets);
	int *tmp2 = (int *) malloc(sizeof(int) * size * buckets);
	int shift = bitsSortedOn, currIt = 0, prevIt, lastUsedTmp = 1;
	// begin by sorting into tmp
	for (long i=0; i < size; i++) {
		long bucket = arr[i] & k;
		tmp1[freq[currIt][bucket] + size * bucket] = arr[i];
		freq[currIt][bucket]++;
	}

	while (shift < 32) {
		prevIt = currIt++;
		// sorts from tmp1 into tmp2
		for (long bucket=0; bucket<buckets; bucket++) {
			long bucketOffset = bucket * size, j = 0;
			long elemsInBucket = freq[prevIt][bucket];
			while (j < elemsInBucket) {
				long elem = tmp1[bucketOffset + j++];
				tmp2[size * ((elem >> shift) & k) + freq[currIt][(elem >> shift) & k]++] = elem;
			}
		}
		shift += bitsSortedOn;
		if (!(shift < 32)) {
			lastUsedTmp = 2;
			break;
		}
		prevIt = currIt++;
		// sorts from tmp2 back into tmp1
		for (long bucket=0; bucket<buckets; bucket++) {
			long bucketOffset = bucket * size, j = 0;
		    long	elemsInBucket = freq[prevIt][bucket];
			while (j < elemsInBucket) {
				long elem = tmp2[bucketOffset + j++];
				tmp1[size * ((elem >> shift) & k) + freq[currIt][(elem >> shift) & k]++] = elem;
			}
		}
		shift += bitsSortedOn;
	}
	shift -= bitsSortedOn;
	int *tmp = (lastUsedTmp == 1) ? tmp1 : tmp2;
	// copy back into array
	for (int bucket=buckets-1, elem = size-1; bucket>=0; bucket--) {
		long bucketOffset = bucket*size;
		while (freq[currIt][bucket] != 0) {
			arr[elem--] = tmp[bucketOffset + --freq[currIt][bucket]];
		}
	}

	free(tmp1);
	free(tmp2);
	return 0;
}

int msdLsdRadixSort(int *arr, int size, int msdBits, int lsdBits) {
	int const msdBuckets = 1 << msdBits;
	int bitMask = msdBuckets-1;
	int *freq = (int *) calloc(msdBuckets, sizeof(int));
	int *tmp = (int *) malloc(sizeof(int) * size);
	// indexes of where each msdBuckets is placed, first one obv begins at 0
	int msdIndexes[msdBuckets+1]; msdIndexes[0] = 0; msdIndexes[msdBuckets] = size;
	int shift = 32 - msdBits;

	// TODO msd-sort
	for (int i=0; i<size; i++) { // count frequencies
		freq[(arr[i] >> shift) & bitMask]++; }
	for (int i=1; i<msdBuckets; i++) { // sum frequencies
		freq[i] += freq[i-1]; 
		msdIndexes[i] = freq[i-1];
	}
	for (int i=size-1; i >= 0; i--) { // move nodes correct loc in tmp array
		int index = --freq[(arr[i] >> shift) & bitMask];
		tmp[index] = arr[i]; 
	}
	for (int i=0; i<size; i++) { // copy from tmp back to arr
		arr[i] = tmp[i]; }

	// variables for sorting lsd sorting
	int const lsdBuckets = 1 << lsdBits;
	bitMask = lsdBuckets-1;
	free(freq); free(tmp);
	freq = (int *) calloc(lsdBuckets, sizeof(int));
	shift = 0;
	int *tmp1 = (int *) malloc(sizeof(int) * size);
	int *tmp2 = (int *) malloc(sizeof(int) * size);

	// TODO sort rest using lsd
	// go through each bucket and sort it in lsd fashion
	for (int bucket=0; bucket < msdBuckets; bucket++) {
		// cout << "sorting bucket " << bucket << endl;
		shift = 0;
		int remainingBits = 32 - msdBits;
		// start and end indexes of the elems in the bucket in arr
		int start = msdIndexes[bucket], end = msdIndexes[bucket+1]; 
		int elemts = end-start;
		bool exitedEarly = false;
		// don't sort if bucket contains less than one element
		if ((end-start) < 2) continue;
		// copy the elems of the bucket into tmp1
		for (int i = start, j = 0; i < end; i++,j++) { 
			tmp1[j] = arr[i];
		}

		// sort the bucket in lsd fashion by sorting between tmp1 and 2
		while (shift < remainingBits) {
			for (int i=0; i<elemts; i++) {
				freq[(tmp1[i] >> shift) & bitMask]++; }
			for (int i=1; i<lsdBuckets; i++) {
				freq[i] += freq[i-1]; }
			for (int i=elemts-1; i >= 0; i--) {
				tmp2[--freq[(tmp1[i] >> shift) & bitMask]] = tmp1[i]; }
			for (int i=0; i<lsdBuckets; i++) {
				freq[i] = 0; }

			shift += lsdBits;
			if (!(shift < remainingBits)) { exitedEarly = true; break; }

			for (int i=0; i<elemts; i++) { // count frequencies
				freq[(tmp2[i] >> shift) & bitMask]++; }
			for (int i=1; i<lsdBuckets; i++) { // sum frequencies
				freq[i] += freq[i-1]; }
			for (int i=elemts-1; i >= 0; i--) { // move nodes to correct loc in tmp1
				tmp1[--freq[(tmp2[i] >> shift) & bitMask]] = tmp2[i]; }
			for (int i=0; i<lsdBuckets; i++) { // set frequencies to 0
				freq[i] = 0; }
			shift += lsdBits;
	   	}
		// copy from the correct bucket depending on which was last sorted into
		tmp = (exitedEarly) ? tmp2 : tmp1;
		// copy back from tmp into the correct position of arr
		for (int i = start, j = 0; i < end; i++,j++) { 
			arr[i] = tmp[j]; }
	}
	free(tmp1); free(tmp2);
	free(freq);
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
	for (int k=0; k<reps; k++) {
		for (int i=1; i<=16;i++) {
			for (int j = 0; j<N; j++) {
				array[j] = rand();
			}
			clock_t t = clock();
			// radixSortWoCountingFreq(array, N, i);
			radixSort(array, N, i);
			times[i] += clock() - t;
		}
	}
	for (int i = 1; i<17; i++) {
		cout << (times[i]/reps) << " when sorting on " << i << " bits" << endl;
	}
	delete[] array;
}

int main() {
}
