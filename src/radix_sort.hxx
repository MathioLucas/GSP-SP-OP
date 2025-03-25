// this file contains a definition of radix sort, used solely in the authentication algorithm of a decomposition tree

#ifndef __RADIX_SORT_HXX__
#define __RADIX_SORT_HXX__

#include <algorithm>

void radix_sort_helper(std::vector<int>& arr, int bit, int start, int end) { // radix sort, which the paper suggests using to sort the adjacency lists of the original graph for comparison with the produced ones from the SP tree in O(|E|) time
										 									 // this implementation is MSD-first, in-place, recursive, with a radix of 16
							 			 									 // I don't particularly like this approach; radix sort doesn't run in O(n) time but O(n * w) time, where w is the maximum length of the keys being sorted
	if (end - start <= 100) { // sort small subarrays with std::sort rather than continuing to recurse (the exact number was decided on a whim)
							  // this is a good idea because std::sort is faster on small arrays due to the lower constant despite being O(nlogn)), while radix sort has quite a hefty constant (we need to initialize all the bucket boundaries, which takes O(1) time but is expensive nonetheless)
		std::sort(arr.begin() + start, arr.begin() + end);
		return;
	}

	int bucket_count[16];
	int bucket_bounds[16];

	for (int i = 0; i < 16; i++) {
		bucket_count[i] = 0;
		bucket_bounds[i] = 0;
	}

	int mask = 15 << bit;
	for (int i = start; i < end; i++) { // count how many of our numbers go into each bucket
		bucket_count[(arr[i] & mask) >> bit]++;
	}

	for (int i = 1; i < 16; i++) bucket_count[i] += bucket_count[i - 1]; // compute prefix sum (so we know where to put our numbers)
	for (int i = 0; i < 16; i++) { // initialize bucket boundaries
		bucket_count[i] += start;
		bucket_bounds[i] = bucket_count[i];
	}

	int i = start;
	for (int curr_gap = 0; curr_gap < 16; curr_gap++) { // partition array into buckets
		while (i < bucket_bounds[curr_gap]) {
			int bucket = (arr[i] & mask) >> bit; // swap this element into its bucket
			bucket_bounds[bucket]--;
			int temp = arr[bucket_bounds[bucket]];
			arr[bucket_bounds[bucket]] = arr[i];
			arr[i] = temp;
		}

		i = bucket_count[curr_gap]; // skip over already partitioned sections
	}

	int substart = start;
	if (bit > 0) {
	 	for (int i = 0; i < 16; i++) { // recurse
			radix_sort_helper(arr, bit - 4, substart, bucket_count[i]);
			substart = bucket_count[i];
		}
	}
}

void radix_sort(std::vector<int>& arr) {
	radix_sort_helper(arr, sizeof(int) * 8 - 4, 0, arr.size());
}

#endif