#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <cstddef>
#include <vector>
#include <sys/stat.h>
#include <errno.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include "variant_permutation_v3.h"

using namespace std;

#define STRSIZE 1000
#define NUMTHREADSBASE 32

// Refactorization of the code that turns a chromosome string into an integer
__host__ int chr2int (string chr_str) {
	if (chr_str == "chrX") {
		return 23;
	} else if (chr_str == "chrY") {
		return 24;
	} else if (chr_str == "chrM" || chr_str == "chrMT") {
		return 25;
	} else {
		string chr_part = chr_str.substr(3);
		return atoi(chr_part.c_str());
	}
}

__device__ bool gpuCmpIntervals (int var_chr, int var_start, int var_end, int ann_chr, int ann_start, int ann_end) {
	if (var_chr != ann_chr) {
		return (var_chr < ann_chr);
	} else if (var_start != ann_start) {
		return (var_start < ann_start);
	} else if (var_end != ann_end) {
		return (var_end < ann_end);
	} else { // The intervals are equal, so return false since a is not less than b
		return false;
	}
}

inline void GPUassert(cudaError_t code, const char * file, int line, bool Abort=true)
{
    if (code != 0) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code),file,line);
        if (Abort) exit(code);
    }       
}

#define GPUerrchk(ans) { GPUassert((ans), __FILE__, __LINE__); }

__device__ void BottomUpMerge(int* target_array, int* temp_array, int left_head, int right_head, int tail) {
	int left_ptr = left_head;
	int right_ptr = right_head;
	int end_ptr = tail;
	
	/* While there are elements in the left or right runs */
  for (int j = left_ptr; j < end_ptr; j++) {
  
  	/* If left run head exists and is <= existing right run head */
  	if (left_ptr < right_head && (right_ptr >= tail || target_array[left_ptr] <= target_array[right_ptr])) {
  		temp_array[j] = target_array[left_ptr];
  		left_ptr++;
  	} else {
  		temp_array[j] = target_array[right_ptr];
  		right_ptr++;
  	}
  }
}

__device__ void CopyArray(int* source, int* dest, int n) {
	for (int i = 0; i < n; i++) {
		dest[i] = source[i];
	}
}

__device__ void BottomUpSort(int* target_array, int n) {
	int *temp_array = (int *)malloc(n*sizeof(int));
	for (int width = 1; width < n; width = 2*width) {
		for (int i = 0; i < n; i = i+2*width) {
			int left_ptr = i;
			int right_ptr = (int)fmin((double)(i+width), (double)n);
			int end_ptr = (int)fmin((double)(i+2*width), (double)n);
			BottomUpMerge(target_array, temp_array, left_ptr, right_ptr, end_ptr);
		}
		// Copy work done in temp_array into target_array for next iteration
		CopyArray(temp_array, target_array, n);
	}
	free(temp_array);
}

__device__ void intersection_kernel(int start, int end, int* gpu_var_chr, int* gpu_var_start, int* gpu_var_end, int* gpu_ann_chr, int* gpu_ann_start, int* gpu_ann_end, int* gpu_var_arr_length, int* gpu_n, int* gpu_dmin, int* gpu_dmax, double* gpu_pvalues) {

	// DEBUG
	// int tid = threadIdx.x + blockIdx.x * blockDim.x;
	// printf("Intersection kernel %d\n", tid);

	for (int i = start; i <= end; i++) {
		// Unpack the current annotation
		int this_ann_chr = gpu_ann_chr[i];
		int this_ann_start = gpu_ann_start[i];
		int this_ann_end = gpu_ann_end[i];
		
		// Initialize the variant pointer
		int vlength = *gpu_var_arr_length;
		vlength = vlength - 1;
		int vthis = vlength/2;
		vlength = vlength/2;
		
		int vlength_const = *gpu_var_arr_length;
		
		int this_var_chr;
		int this_var_start;
		int this_var_end;
		
		// Keep track of whether the target is greater or less than the current variant
		// Also what was the last comparison
		int prev_is_greater = -1;
		int is_greater = -1;
		
		int int_variants;
		
		// return;
		// DEBUG
		// int test = 0;
		
		while (1) {
			// Unpack current variant
			this_var_chr = gpu_var_chr[vthis];
			this_var_start = gpu_var_start[vthis];
			this_var_end = gpu_var_end[vthis];
			
			// DEBUG
// 			printf("Variant %d: %d, %d, %d\n", test, this_var_chr, this_var_start, this_var_end);
// 			test++;
		
			// Check for intersection
			if (this_var_chr == this_ann_chr && this_var_start <= this_ann_end && this_ann_start <= this_var_end) {
				int_variants = 1;
				break;
			} else {
				if (vlength > 1) { // vlength does not fall below 1
					vlength = vlength/2;
				}
				if (!(gpuCmpIntervals(this_var_chr, this_var_start, this_var_end, this_ann_chr, this_ann_start, this_ann_end))) {
					if (vthis == 0) { // Don't read off the end of the array
						int_variants = 0;
						break;
					} else {
						// Take the smaller half
						vthis = vthis-vlength;
						prev_is_greater = is_greater;
						is_greater = 0;
					}
				} else {
					if (vthis == vlength_const - 1) { // Don't read off the end of the array
						int_variants = 0;
						break;
					} else {
						// Take the larger half
						vthis = vthis+vlength;
						prev_is_greater = is_greater;
						is_greater = 1;
					}
				}
				if (vlength == 1 && ((prev_is_greater == 1 && is_greater == 0) || (prev_is_greater == 0 && is_greater == 1))) { // No intersection
					int_variants = 0;
					break;
				}
			}
		}
		
		// DEBUG
		// printf("After intersection found\n");
		// printf("<-- %d, %d, %d -->\n", this_ann_chr, this_ann_start, this_ann_end);
		// printf("int_variants: %d\n", int_variants);
		
		// int int_variants = 1;
		
		int v_anchor = vthis;
		
		if (v_anchor != 0) {
			vthis--;
		
			// Unpack current variant
			this_var_chr = gpu_var_chr[vthis];
			this_var_start = gpu_var_start[vthis];
			this_var_end = gpu_var_end[vthis];
		
			// Search for intersecting variants bidirectionally
			while (this_var_chr == this_ann_chr && this_var_start <= this_ann_end && this_ann_start <= this_var_end) {
				int_variants++;
				vthis--;
				this_var_chr = gpu_var_chr[vthis];
				this_var_start = gpu_var_start[vthis];
				this_var_end = gpu_var_end[vthis];
			}
		}
		
		if (v_anchor != vlength_const-1) {
			vthis = v_anchor;
			vthis++;
		
			// Unpack current variant
			this_var_chr = gpu_var_chr[vthis];
			this_var_start = gpu_var_start[vthis];
			this_var_end = gpu_var_end[vthis];
		
			// Search for intersecting variants bidirectionally
			while (this_var_chr == this_ann_chr && this_var_start <= this_ann_end && this_ann_start <= this_var_end) {
				int_variants++;
				vthis++;
				this_var_chr = gpu_var_chr[vthis];
				this_var_start = gpu_var_start[vthis];
				this_var_end = gpu_var_end[vthis];
			}
		}
		
		// DEBUG
		// printf("Find random bins\n");
		// printf("<-- %d, %d, %d -->\n", this_ann_chr, this_ann_start, this_ann_end);
		// printf("int_variants: %d\n", int_variants);
		
		// Number of random bins to select
		int n = (*gpu_n);
	
		// The minimum distance between element and random bin
		int dmin = (*gpu_dmin);
	
		// The maximum distance between element and random bin
		int dmax = (*gpu_dmax);
		
		// Pick random bins from surrounding regions
		// We take n/2 from the upstream region, and n/2 from the downstream regions
		// Random number drawn from [0, dmax - dmin - (annotation's length)]
		// Save only the start coordinate, the chr and end can be derived JIT
		int range = dmax - dmin - (this_ann_end - this_ann_start + 1);
		int *upstream_start = (int *)malloc((n/2)*sizeof(int));
		int *downstream_start = (int *)malloc((n/2)*sizeof(int));
		
		// Upstream bin selection
		// Configure where the start of this range is
		// int rand_range_chr = this_ann_chr;
		int rand_range_start = this_ann_start - dmax;
		
		curandState *d_state;
		d_state = (curandState *)malloc(sizeof(curandState));
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		curand_init(65536, tid, 0, d_state);
		
		for (int j = 0; j < n/2; j++) {
			float this_rand = curand_uniform(d_state);
			int rand_start = this_rand*range;
			rand_start += rand_range_start;
			// int rand_start = rand() % range;
			
			// DEBUG
			// printf("this_rand_start: %d\n", rand_start);
			
			upstream_start[j] = rand_start;
		}
		
		// Downstream bin selection
		// Configure where the start of this range is
		// rand_range_chr = this_ann_chr;
		rand_range_start = this_ann_end + dmin;
		
		for (int j = 0; j < n/2; j++) {
			float this_rand = curand_uniform(d_state);
			int rand_start = this_rand*range;
			rand_start += rand_range_start;
			// int rand_start = rand() % range;
			
			// DEBUG
			// printf("this_rand_start: %d\n", rand_start);
			downstream_start[j] = rand_start;
		}
		
		// DEBUG: Check upstream and downstream random bins
// 		for (int k = 0; k < (n/2); k++) {
// 			printf("Upstream bin %d: %d\n", k, upstream_start[k]);
// 		}
// 		for (int k = 0; k < (n/2); k++) {
// 			printf("Downstream bin %d: %d\n", k, downstream_start[k]);
// 		}
		
		// Find the intersecting variants for the random bins
		
		// Sort the upstream and downstream random bins
		BottomUpSort(upstream_start, (n/2));
		BottomUpSort(downstream_start, (n/2));
		// gpuIntervalMergeSortByEnd(upstream_chr, upstream_start, upstream_end, upstream_chr_sorted, upstream_start_sorted, upstream_end_sorted, 0, (n/2)-1);
		// gpuIntervalMergeSortByStart(downstream_chr, downstream_start, downstream_end, downstream_chr_sorted, downstream_start_sorted, downstream_end_sorted, 0, (n/2)-1);
		
		// DEBUG: Check upstream and downstream random bins (sorted)
// 		for (int k = 0; k < (n/2); k++) {
// 			printf("Upstream bin (sorted) %d: %d\n", k, upstream_start[k]);
// 		}
// 		for (int k = 0; k < (n/2); k++) {
// 			printf("Downstream bin (sorted) %d: %d\n", k, downstream_start[k]);
// 		}
		
	 	// Upstream bins: search backwards from the variant at gpu_var_*[v_anchor]
	 	unsigned int vpointer2 = v_anchor;
	 	
	 	// A collection of intersecting variants counts from the random bins
	 	int *varcounts = (int *)malloc(n*sizeof(int));
	 	int varcounts_index = 0;
	 	
	 	// Backwards search!
	 	unsigned int j = (n/2);
	 	do {
	 		j--;
	 		
	 		// How many variants intersect this bin?
 			int this_variants = 0;
 			
 			// Unpack the current annotation
 			int upstream_ann_chr = this_ann_chr;
 			int upstream_ann_start = upstream_start[j];
 			int upstream_ann_end = upstream_ann_start + (this_ann_end - this_ann_start);
 			
 			// Unpack the current variant
 			int upstream_var_chr = gpu_var_chr[vpointer2];
 			int upstream_var_start = gpu_var_start[vpointer2];
 			int upstream_var_end = gpu_var_end[vpointer2];
 			
 			// Now count the intersecting variants
 			// vpointer2 points to the "earliest" possible annotation, and vpointer3
			// points to the variants up until the last intersecting with the annotation
			unsigned int vpointer3 = vpointer2;
			
			// While vpointer3 does not go past the current annotation
			while (upstream_var_chr > upstream_ann_chr || (upstream_var_chr == upstream_ann_chr && upstream_var_start >= upstream_ann_start)) {
				
				// If the current variant intersects the current annotation, increment target_variants
				if (upstream_var_chr == upstream_ann_chr && upstream_ann_start <= upstream_var_end && upstream_var_start <= upstream_ann_end) {
					this_variants++;
				} else { // Update vpointer2
					if (vpointer3 != 0) {
						vpointer2 = vpointer3 - 1;
					}
				}
				// Now update the cur_var
				if (vpointer3 == 0) {
					break;
				}
				vpointer3--;
				
				upstream_var_chr = gpu_var_chr[vpointer3];
				upstream_var_start = gpu_var_start[vpointer3];
				upstream_var_end = gpu_var_end[vpointer3];
			}
			
			// this_variants has been settled, save for output
			varcounts[varcounts_index] = this_variants;
			varcounts_index++;
		} while (j > 0);
		
		// Downstream bins: a more straight forward search :)
		vpointer2 = v_anchor;
		
		for (unsigned int j = 0; j < (n/2); j++) {
			
			// How many variants intersect this bin?
			int this_variants = 0;
			
			// Unpack the current annotation
			int downstream_ann_chr = this_ann_chr;
 			int downstream_ann_start = downstream_start[j];
 			int downstream_ann_end = downstream_ann_start + (this_ann_end - this_ann_start);
 			
 			// Unpack the current variant
 			int downstream_var_chr = gpu_var_chr[vpointer2];
 			int downstream_var_start = gpu_var_start[vpointer2];
 			int downstream_var_end = gpu_var_end[vpointer2];
 			
 			// Now count the intersecting variants
 			// vpointer2 points to the "earliest" possible variant, and vpointer3
 			// points to the variants up until the last intersecting with the annotation
 			unsigned int vpointer3 = vpointer2;
 			
 			// While vpointer3 does not go past the current annotation
 			while (downstream_var_chr < downstream_ann_chr || (downstream_var_chr == downstream_ann_chr && downstream_var_end <= downstream_ann_end)) {
 				
 				// If the current variant intersects the current annotation, increment target_variants
 				if (downstream_var_chr == downstream_ann_chr && downstream_ann_start <= downstream_var_end && downstream_var_start <= downstream_ann_end) {
 					this_variants++;
 				} else { // Update vpointer2
 					if (vpointer3 != (vlength_const)-1) {
 						vpointer2 = vpointer3 + 1;
 					}
 				}
 				// Now update the cur_var
 				if (vpointer3 == (vlength_const)-1) {
 					break;
 				}
 				vpointer3++;
 				
 				downstream_var_chr = gpu_var_chr[vpointer3];
 				downstream_var_start = gpu_var_start[vpointer3];
 				downstream_var_end = gpu_var_end[vpointer3];
 			}
 			
 			// this_variants has been settled, save for output
 			varcounts[varcounts_index] = this_variants;
 			varcounts_index++;
 		}
 		
 		// DEBUG
//  		for (int k = 0; k < n; k++) {
//  			printf("Varcounts %d: %d\n", k, varcounts[k]);
//  		}
 		
 		// P-value calculation: how many of the random bins have at least as many
 		// variants at k_t?
 		int overbins = 0;
 		for (unsigned int j = 0; j < n; j++) {
 			if (varcounts[j] >= int_variants) {
 				overbins++;
 			}
 		}
 		
 		double fraction = (double)overbins/(double)n;
 		gpu_pvalues[i] = fraction;
 		
 		// Malloc free the temp arrays
 		free(upstream_start);
 		free(downstream_start);
 		free(varcounts);
 		free(d_state);
 		
 		// DEBUG
 		// printf("GPU pvalue %d: %f\n", i, fraction);
	}
}

__global__ void apportionWork(int* gpu_var_chr, int* gpu_var_start, int* gpu_var_end, int* gpu_ann_chr, int* gpu_ann_start, int* gpu_ann_end, int* gpu_var_arr_length, int* gpu_ann_arr_length, int* gpu_n, int* gpu_dmin, int* gpu_dmax, double *gpu_pvalues) {
// __global__ void apportionWork() {

	// DEBUG
	// printf("Running thread\n");
	// *test_int_gpu = 247;
	
	// Which thread am I?
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int total_threads = NUMTHREADSBASE*NUMTHREADSBASE;
	
	// DEBUG
	// printf("%d\n", *gpu_ann_arr_length);
	// printf("%d\n", total_threads);
	
	int length = *gpu_ann_arr_length;
	
	// int kernel_annsize = (*gpu_ann_arr_length)/4;
	int kernel_annsize = length/total_threads;
	int mod = length%total_threads;
	int start;
	int end;
	
	if (kernel_annsize > 0) {
		if (tid < mod) {
			start = tid*(kernel_annsize+1);
			end = ((tid+1)*(kernel_annsize+1))-1;
		} else {
			start = (mod*(kernel_annsize+1))+((tid-mod)*kernel_annsize);
			end = ((mod*(kernel_annsize+1))+((tid-mod+1)*kernel_annsize))-1;
		}
		
		// DEBUG: Print the thread ID, start index, and end index
// 		printf("Thread ID: %d; start index: %d; end index: %d\n", tid, start, end);
// 		return;
		
		intersection_kernel(start, end, gpu_var_chr, gpu_var_start, gpu_var_end, gpu_ann_chr, gpu_ann_start, gpu_ann_end, gpu_var_arr_length, gpu_n, gpu_dmin, gpu_dmax, gpu_pvalues);
	} else {
		start = tid;
		end = tid;
		
		// DEBUG: Print the thread ID, start index, and end index
		// printf("Thread ID: %d; start index: %d; end index: %d\n", tid, start, end);
		
		if (tid < length) {
			intersection_kernel(start, end, gpu_var_chr, gpu_var_start, gpu_var_end, gpu_ann_chr, gpu_ann_start, gpu_ann_end, gpu_var_arr_length, gpu_n, gpu_dmin, gpu_dmax, gpu_pvalues);
		}
	}
}

/*
 * Subroutine that merges the input intervals (3col)
 * Assumes input interval array is sorted
 */
__host__ vector<vector<string> > merge_intervals (vector<vector<string> > starting_intervals) {
	
	// Output vector
	vector<vector<string> > resulting_intervals;
	
	for (unsigned int i = 0; i < starting_intervals.size(); i++) {
		if (i == 0) { // Just push the first interval straight onto resulting_intervals
			resulting_intervals.push_back(starting_intervals[i]);
		} else { // Compare the i'th starting to the last resulting
			vector<string> interval1 = resulting_intervals[resulting_intervals.size()-1];
			vector<string> interval2 = starting_intervals[i];
			
			vector<vector<string> > vec1;
			vec1.push_back(interval1);
			
			vector<vector<string> > int_interval = intersecting_intervals(vec1, interval2);
			
			// If there's anything in int_interval, there is an overlap, otherwise no
			if (int_interval.size() == 0) { // No intersection
				resulting_intervals.push_back(starting_intervals[i]);
			} else { // Yes intersection, create a merged interval
				vector<string> merged;
				merged.push_back(interval1[0]);
				int new_start = min(atoi(interval1[1].c_str()), atoi(interval2[1].c_str()));
				int new_end = max(atoi(interval1[2].c_str()), atoi(interval2[2].c_str()));
				
				char new_start_cstr[STRSIZE];
				sprintf(new_start_cstr, "%d", new_start);
				merged.push_back(string(new_start_cstr));
		
				char new_end_cstr[STRSIZE];
				sprintf(new_end_cstr, "%d", new_end);
				merged.push_back(string(new_end_cstr));
				
				resulting_intervals.pop_back();
				resulting_intervals.push_back(merged);
			}
		}
	}
	return resulting_intervals;
}

/* 
 * This code takes as input a variant track and an annotation track. For each
 * element, count intersecting variants, and select n (hardcoded) random bins
 * from a range upstream and downstream of the element. These random bins will
 * fall within d_min and d_max distance from the element, and the n bins will
 * be split evenly between upstream and downstream bins. Then, calculate a p-value
 * for mutation burden based on how many of these random bins have more intersecting
 * variants.
 */
 
int main (int argc, char* argv[]) {
	
	/* User-supplied arguments */
	
	// Number of random bins to select
	int n;
	
	// The minimum distance between element and random bin
	int dmin;
	
	// The maximum distance between element and random bin
	int dmax;
	
	// File with prohibited coordinates
	// Expected format: tab(chr, start, end, ...)
	string prohibited_file;
	
	// File with single nucleotide variants
	// Expected format: tab(chr, start, end, ...)
	string vfile;
	
	// File with annotations to study for mutation burden
	// Expected format: tab(chr, start, end, name, ...)
	string afile;
	
	// File with the output
	// Format: tab(chr, start, end, name, p-value)
	string outfile;
	
	if (argc != 8) {
		fprintf(stderr, "Usage: moat_a_gpu [# of permutations] [d_min] [d_max] [prohibited regions file] [variant file] [annotation file] [output file]. Exiting.\n");
		return 1;
	} else {
		n = atoi(argv[1]);
		dmin = atoi(argv[2]);
		dmax = atoi(argv[3]);
		prohibited_file = string(argv[4]);
		vfile = string(argv[5]);
		afile = string(argv[6]);
		outfile = string(argv[7]);
	}
	
	// Verify files, and import data to memory
	struct stat vbuf;
	if (stat(vfile.c_str(), &vbuf)) { // Report the error and exit
		fprintf(stderr, "Error trying to stat %s: %s\n", vfile.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (vbuf.st_size == 0) {
		fprintf(stderr, "Error: Variant file cannot be empty. Exiting.\n");
		return 1;
	}
	
	struct stat abuf;
	if (stat(afile.c_str(), &abuf)) { // Report the error and exit
		fprintf(stderr, "Error trying to stat %s: %s\n", afile.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (abuf.st_size == 0) {
		fprintf(stderr, "Error: Annotation file cannot be empty. Exiting.\n");
		return 1;
	}
	
	struct stat pbuf;
	if (stat(prohibited_file.c_str(), &pbuf)) { // Report the error and exit
		fprintf(stderr, "Error trying to stat %s: %s\n", prohibited_file.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (pbuf.st_size == 0) {
		fprintf(stderr, "Error: Prohibited regions file cannot be empty. Exiting.\n");
		return 1;
	}
	
	/* Data structures for the starting data */
	// Variant arrays, contains variants of the format vector(chr, start, end)
	vector<vector<string> > var_array;
	
	// Annotation arrays, contains annotations of the format vector(chr, start, end, name)
	vector<vector<string> > ann_array;
	
	// Prohibited regions array, contains annotations of the format vector(chr, start, end)
	vector<vector<string> > prohibited_regions;
	
	// DEBUG
	// printf("Breakpoint 1\n");
	
	// Bring variant file data into memory
	// Save the first 3 columns, ignore the rest if there are any
	char linebuf[STRSIZE];
	FILE *vfile_ptr = fopen(vfile.c_str(), "r");
	while (fgets(linebuf, STRSIZE, vfile_ptr) != NULL) {
		string line = string(linebuf);
		
		// DEBUG
		// printf("%s\n", line.c_str());
		
		// Extract chromosome, start, and end from line (first 3 columns)
		vector<string> vec;
		for (int i = 0; i < 3; i++) {
			size_t ws_index = line.find_first_of("\t\n");
			string in = line.substr(0, ws_index);
			vec.push_back(in);
			line = line.substr(ws_index+1);
		}
		
		// If this is not a standard chromosome, then remove this row
		if (chr2int(vec[0]) == 0) {
			continue;
		}
		
		var_array.push_back(vec);
	}
	// Check feof of vfile
	if (feof(vfile_ptr)) { // We're good
		fclose(vfile_ptr);
	} else { // It's an error
		char errstring[STRSIZE];
		sprintf(errstring, "Error reading from %s", vfile.c_str());
		perror(errstring);
		return 1;
	}
	
	// Bring annotation file data into memory
	FILE *afile_ptr = fopen(afile.c_str(), "r");
	while (fgets(linebuf, STRSIZE, afile_ptr) != NULL) {
		string line = string(linebuf);
		
		// DEBUG
		// printf("%s", line.c_str());
		
		// Extract chromosome, start, end, and name from line (first 4 columns)
		vector<string> vec;
		for (int i = 0; i < 4; i++) {
			size_t ws_index = line.find_first_of("\t\n");
			string in = line.substr(0, ws_index);
			vec.push_back(in);
			line = line.substr(ws_index+1);
		}
		
		// If this is not a standard chromosome, then remove this row
		if (chr2int(vec[0]) == 0) {
			continue;
		}
		
		ann_array.push_back(vec);
	}
	// Check feof of vfile
	if (feof(afile_ptr)) { // We're good
		fclose(afile_ptr);
	} else { // It's an error
		char errstring[STRSIZE];
		sprintf(errstring, "Error reading from %s", afile.c_str());
		perror(errstring);
		return 1;
	}
	
	// Import prohibited regions file
	FILE *prohibited_file_ptr = fopen(prohibited_file.c_str(), "r");
	while (fgets(linebuf, STRSIZE, prohibited_file_ptr) != NULL) {
	
		string line = string(linebuf);
		
		// Extract chromosome, start, and end from line (first 3 columns)
		vector<string> vec;
		for (int i = 0; i < 3; i++) {
			size_t ws_index = line.find_first_of("\t\n");
			string in = line.substr(0, ws_index);
			vec.push_back(in);
			line = line.substr(ws_index+1);
		}
		
		// If this is not a standard chromosome, then remove this row
		if (chr2int(vec[0]) == 0) {
			continue;
		}
		
		prohibited_regions.push_back(vec);
	}
	// Check feof of prohibited_file_ptr
	if (feof(prohibited_file_ptr)) { // We're good
		fclose(prohibited_file_ptr);
	} else { // It's an error
		char errstring[STRSIZE];
		sprintf(errstring, "Error reading from %s", prohibited_file.c_str());
		perror(errstring);
		return 1;
	}
	
	// DEBUG
	// printf("Breakpoint 2\n");
	
	// Sort the arrays
	sort(var_array.begin(), var_array.end(), cmpIntervals);
	sort(ann_array.begin(), ann_array.end(), cmpIntervals);
	sort(prohibited_regions.begin(), prohibited_regions.end(), cmpIntervals);
	
	// Merge prohibited regions
	prohibited_regions = merge_intervals(prohibited_regions);
	
	// Remove variants and annotations that intersect the blacklist regions
// 	vector<vector<string> > var_array_new;
	for (unsigned int i = 0; i < var_array.size(); i++) {
		vector<vector<string> > inter = intersecting_intervals(prohibited_regions, var_array[i]);
		if (inter.size() > 0) {
			// var_array_new.push_back(var_array[i]);
			var_array[i][0] = "chrNo";
		}
	}
// 	var_array = var_array_new;

	sort(var_array.begin(), var_array.end(), cmpIntervals);
	
	// Remove those marked for deletion
	while (var_array[var_array.size()-1][0] == "chrNo") {
		var_array.erase(var_array.end());
	}
	
	// vector<vector<string> > ann_array_new;
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		vector<vector<string> > inter = intersecting_intervals(prohibited_regions, ann_array[i]);
		if (inter.size() > 0) {
			// ann_array_new.push_back(ann_array[i]);
			ann_array[i][0] = "chrNo";
		}
	}
	// ann_array = ann_array_new;
	
	sort(ann_array.begin(), ann_array.end(), cmpIntervals);
	
	// Remove those marked for deletion
	while (ann_array[ann_array.size()-1][0] == "chrNo") {
		ann_array.erase(ann_array.end());
	}
	
	// DEBUG
// 	printf("DEBUG: Sorted var_array\n");
// 	for (unsigned int i = 0; i < var_array.size(); i++) {
// 		printf("%s, %s, %s\n", var_array[i][0].c_str(), var_array[i][1].c_str(), var_array[i][2].c_str());
// 	}
	// printf("Breakpoint 3\n");
	
	// Variables for main loop
	// unsigned int var_pointer = 0;
	
	// Can malloc free the prohibited regions
	prohibited_regions.clear();
	
	// Length of each variant array
	int var_arr_length = var_array.size();
	
	// Turn the vectors into int arrays for CUDA
	// Var array
	int *var_chr = (int*)malloc(var_arr_length*sizeof(int));
	int *var_start = (int*)malloc(var_arr_length*sizeof(int));
	int *var_end = (int*)malloc(var_arr_length*sizeof(int));
	
	// Lengths of each annotation array
	int ann_arr_length = ann_array.size();
	
	// Ann array
	int *ann_chr = (int*)malloc(ann_arr_length*sizeof(int));
	int *ann_start = (int*)malloc(ann_arr_length*sizeof(int));
	int *ann_end = (int*)malloc(ann_arr_length*sizeof(int));
	
	// Variant array processing
	for (unsigned int i = 0; i < var_array.size(); i++) {
		// Unpack the current variant
		string cur_var_chr = var_array[i][0];
		string cur_var_start = var_array[i][1];
		string cur_var_end = var_array[i][2];
		
		int cur_var_chr_int;
		if (cur_var_chr == "chrX") {
			cur_var_chr_int = 24;
		} else if (cur_var_chr == "chrY") {
			cur_var_chr_int = 25;
		} else if (cur_var_chr == "chrM") {
			cur_var_chr_int = 26;
		} else {
			string cur_var_chr_part = cur_var_chr.substr(3);
			cur_var_chr_int = atoi(cur_var_chr_part.c_str());
		}
		
		// var_arr_length++;
		
		// var_chr = (int*)realloc(var_chr, var_arr_length*sizeof(int));
		var_chr[i] = cur_var_chr_int;
		
		// var_start = (int*)realloc(var_start, var_arr_length*sizeof(int));
		var_start[i] = atoi(cur_var_start.c_str());
		
		// var_end = (int*)realloc(var_end, var_arr_length*sizeof(int));
		var_end[i] = atoi(cur_var_end.c_str());
	}
	
	// Can malloc free the variant array
	var_array.clear();
	
	// Annotation array processing
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		// Unpack the current variant
		string cur_ann_chr = ann_array[i][0];
		string cur_ann_start = ann_array[i][1];
		string cur_ann_end = ann_array[i][2];
		
		int cur_ann_chr_int;
		if (cur_ann_chr == "chrX") {
			cur_ann_chr_int = 24;
		} else if (cur_ann_chr == "chrY") {
			cur_ann_chr_int = 25;
		} else if (cur_ann_chr == "chrM") {
			cur_ann_chr_int = 26;
		} else {
			string cur_ann_chr_part = cur_ann_chr.substr(3);
			cur_ann_chr_int = atoi(cur_ann_chr_part.c_str());
		}
		
		// ann_arr_length++;
		
		// ann_chr = (int*)realloc(ann_chr, ann_arr_length*sizeof(int));
		ann_chr[i] = cur_ann_chr_int;
		
		// ann_start = (int*)realloc(ann_start, ann_arr_length*sizeof(int));
		ann_start[i] = atoi(cur_ann_start.c_str());
		
		// ann_end = (int*)realloc(ann_end, ann_arr_length*sizeof(int));
		ann_end[i] = atoi(cur_ann_end.c_str());
	}
	
	// Can't malloc free the annotation array
	// ann_array.clear();
	
	// DEBUG
	// printf("Breakpoint 4\n");
	
	// Begin the CUDA magic
	int *gpu_var_chr;
	int *gpu_var_start;
	int *gpu_var_end;
	
	int *gpu_ann_chr;
	int *gpu_ann_start;
	int *gpu_ann_end;
	
	int *gpu_var_arr_length;
	int *gpu_ann_arr_length;
	
	int *gpu_n;
	int *gpu_dmin;
	int *gpu_dmax;
	
	// DEBUG
// 	printf("Begin CUDA code\n");
// 	int *test_int_cpu;
// 	int *test_int_gpu;
// 	test_int_cpu = (int*)malloc(sizeof(int));
// 	cudaMalloc((void**)&test_int_gpu, sizeof(int));
// 	*test_int_cpu = 246;
// 	cudaMemcpy(test_int_gpu, test_int_cpu, sizeof(int), cudaMemcpyHostToDevice);
	
	cudaMalloc((void**)&gpu_var_chr, var_arr_length*sizeof(int));
	// GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_var_start, var_arr_length*sizeof(int));
	// GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_var_end, var_arr_length*sizeof(int));
	// GPUerrchk(cudaPeekAtLastError());
	
	cudaMalloc((void**)&gpu_ann_chr, ann_arr_length*sizeof(int));
	// GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_ann_start, ann_arr_length*sizeof(int));
	// GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_ann_end, ann_arr_length*sizeof(int));
	// GPUerrchk(cudaPeekAtLastError());
	
	cudaMalloc((void**)&gpu_var_arr_length, sizeof(int));
	// GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_ann_arr_length, sizeof(int));
	// GPUerrchk(cudaPeekAtLastError());
	
	cudaMalloc((void**)&gpu_n, sizeof(int));
	cudaMalloc((void**)&gpu_dmin, sizeof(int));
	cudaMalloc((void**)&gpu_dmax, sizeof(int));
	
	double *gpu_pvalues;
	cudaMalloc((void**)&gpu_pvalues, ann_arr_length*sizeof(double));
	// GPUerrchk(cudaPeekAtLastError());
	
	cudaMemcpy(gpu_var_chr, var_chr, var_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	// GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_var_start, var_start, var_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	// GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_var_end, var_end, var_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	// GPUerrchk(cudaPeekAtLastError());

	cudaMemcpy(gpu_ann_chr, ann_chr, ann_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	// GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_ann_start, ann_start, ann_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	// GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_ann_end, ann_end, ann_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	// GPUerrchk(cudaPeekAtLastError());
	
	cudaMemcpy(gpu_var_arr_length, &var_arr_length, sizeof(int), cudaMemcpyHostToDevice);
	// GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_ann_arr_length, &ann_arr_length, sizeof(int), cudaMemcpyHostToDevice);
	// GPUerrchk(cudaPeekAtLastError());
	
	cudaMemcpy(gpu_n, &n, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_dmin, &dmin, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_dmax, &dmax, sizeof(int), cudaMemcpyHostToDevice);
	
	// Try out 16x16 and see how that goes
	int num_blocks = NUMTHREADSBASE;
	int threads_per_block = NUMTHREADSBASE;
	
	// DEBUG
	// printf("This is debug print 1\n");
	
	// DEBUG
	// printf("Breakpoint 5\n");
	
	// Adjust the heap size based on the size of the dataset
// 	if (ann_array.size() > 3000) {
// 		int fold = ann_array.size()/3000;
// 		int new_heap_size = fold*8000000;
// 		int new_heap_size = 800000000;
//  		cudaDeviceSetLimit(cudaLimitMallocHeapSize, new_heap_size);
// 	}
	
	apportionWork<<<num_blocks, threads_per_block>>>(gpu_var_chr, gpu_var_start, gpu_var_end, gpu_ann_chr, gpu_ann_start, gpu_ann_end, gpu_var_arr_length, gpu_ann_arr_length, gpu_n, gpu_dmin, gpu_dmax, gpu_pvalues);
	GPUerrchk(cudaPeekAtLastError());
	// apportionWork<<<1,1>>>();
	
	// DEBUG
// 	printf("This is debug print 2\n");
// 	cudaMemcpy(test_int_cpu, test_int_gpu, sizeof(int), cudaMemcpyDeviceToHost);
// 	printf("Test int: %d\n", *test_int_cpu);

	// DEBUG
	// printf("Breakpoint 6\n");
	
	// GPUerrchk(cudaDeviceSynchronize());
	
	// Collect the output values, will end with same size as ann_array
	double *pvalues = (double *)malloc(ann_array.size()*sizeof(double));
	// int block = 1000;
	
	cudaMemcpy(pvalues, gpu_pvalues, ann_arr_length*sizeof(double), cudaMemcpyDeviceToHost);
	// GPUerrchk(cudaPeekAtLastError());
//  	if (gpu_pvalues == NULL) {
//  		printf("Malloc fail!\n");
//  	}
// 	return 0;
// 	
// 	if (ann_arr_length < block) {
// 		cudaMemcpy(pvalues, gpu_pvalues, ann_arr_length*sizeof(double), cudaMemcpyDeviceToHost);
// 		GPUerrchk(cudaPeekAtLastError());
// 	} else {
// 		double *pvalues_ptr = pvalues;
// 		double *gpu_pvalues_ptr = gpu_pvalues;
// 		for (int k = 0; k < ann_arr_length; k += block) {
// 		
// 			// DEBUG
// 			printf("k: %d\n", k);
// 		
// 			int copyblock;
// 			if (k < ann_arr_length-block) {
// 				copyblock = ann_arr_length-k;
// 			} else {
// 				copyblock = block;
// 			}
// 			cudaMemcpy(pvalues_ptr, gpu_pvalues_ptr, copyblock*sizeof(double), cudaMemcpyDeviceToHost);
// 			GPUerrchk(cudaPeekAtLastError());
// 			pvalues_ptr += block*sizeof(double);
// 			gpu_pvalues_ptr += block*sizeof(double);
// 		}
// 	}
	// GPUerrchk(cudaDeviceSynchronize());
	
	// DEBUG
	// printf("Passed CUDA memcpy\n");
// 	for (int i = 0; i < ann_arr_length; i++) {
// 		printf("Pvalue %d: %f\n", i, gpu_pvalues[i]);
// 	}

	// DEBUG
	// printf("Breakpoint 7\n");

	// Output generation
	FILE *outfile_ptr = fopen(outfile.c_str(), "w");
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		
		// Unpack the annotation
		string cur_ann_chr = ann_array[i][0];
		string cur_ann_start = ann_array[i][1];
		string cur_ann_end = ann_array[i][2];
		string cur_ann_name = ann_array[i][3];
		
		// Print the output line
		fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%f\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str(), pvalues[i]);
	}
	fclose(outfile_ptr);
	
 	return 0;
}
