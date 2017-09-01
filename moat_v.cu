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
#include <stdexcept>
#include "variant_permutation_v3.h"

using namespace std;

#define STRSIZE 1000
// #define NUMTHREADSBASE 70
// #define NUMTHREADSBASE 55
// #define NUMTHREADSBASE 32
// #define NUMTHREADSBASE 21
// #define NUMTHREADSBASE 15

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

// Reverse the direction of chr2int
__host__ string int2chr (int chr_num) {
	if (chr_num == 23) {
		return "chrX";
	} else if (chr_num == 24) {
		return "chrY";
	} else if (chr_num == 25) {
		return "chrM";
	} else {
		char chr[STRSIZE];
		sprintf(chr, "%d", chr_num);
		return "chr" + string(chr);
	}
}

__device__ bool gpuCmpIntervals (int var_chr, int var_end, int ann_chr, int ann_start, int ann_end) {
	if (var_chr != ann_chr) {
		return (var_chr < ann_chr);
	} else {
		return (var_end <= ann_end);
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

__device__ void BottomUpSort(int* target_array, int n, int* temp_array) {
	// int *temp_array = (int *)malloc(n*sizeof(int));
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
	// free(temp_array);
}

__device__ void intersection_kernel(int start, int end, int* gpu_var_chr, 
	int* gpu_var_end, double* gpu_var_signal, int* gpu_ann_chr, int* gpu_ann_start, 
	int* gpu_ann_end, int* gpu_var_arr_length, int* gpu_trimer, 
	double* gpu_permuted_chr, double *gpu_permuted_end, int *gpu_wg_switch, 
	curandState **d_state, char *gpu_fasta_1, char *gpu_fasta_2, 
	char *gpu_fasta_3, char *gpu_fasta_4, char *gpu_fasta_5, char *gpu_fasta_6, 
	char *gpu_fasta_7, char *gpu_fasta_8, char *gpu_fasta_9, char *gpu_fasta_10, 
	char *gpu_fasta_11, char *gpu_fasta_12, char *gpu_fasta_13, char *gpu_fasta_14, 
	char *gpu_fasta_15, char *gpu_fasta_16, char *gpu_fasta_17, char *gpu_fasta_18, 
	char *gpu_fasta_19, char *gpu_fasta_20, char *gpu_fasta_21, char *gpu_fasta_22, 
	char *gpu_fasta_23, char *gpu_fasta_24, char *gpu_fasta_25, int *pos_array) {

	// DEBUG
	// int tid = threadIdx.x + blockIdx.x * blockDim.x;
	// printf("Intersection kernel %d\n", tid);
	// printf("GPU var signal 0: %f\n", gpu_var_signal[0]);
	
	// Hardcoded chromosome lengths
	int hg19_coor[25] = {249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 
											 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 
											 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 
											 59128983, 63025520, 48129895, 51304566, 155270560, 59373566, 16571};
	
	int funseq_opt = (*gpu_wg_switch);

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
		// int this_var_start;
		int this_var_end;
		double this_var_signal;
		
		// Keep track of whether the target is greater or less than the current variant
		// Also what was the last comparison
		int prev_is_greater = -1;
		int is_greater = -1;
		
		// int int_variants;
		// Left and right pointers indicate the range of variants that intersect this
		// annotation
		int left_pointer = -1;
		int right_pointer = -1;
		double signal_sum = 0.0;
		
		// return;
		// DEBUG
		// int test = 0;
		
		while (1) {
			// Unpack current variant
			this_var_chr = gpu_var_chr[vthis];
			// this_var_start = gpu_var_start[vthis];
			this_var_end = gpu_var_end[vthis];
			
			if (funseq_opt) {
				this_var_signal = gpu_var_signal[vthis];
			}
			
			// DEBUG
// 			printf("Variant %d: %d, %d, %d\n", test, this_var_chr, this_var_start, this_var_end);
// 			test++;
		
			// Check for intersection
			if (this_var_chr == this_ann_chr && this_ann_start < this_var_end && this_var_end <= this_ann_end) {
				// int_variants = 1;
				left_pointer = vthis;
				right_pointer = vthis;
				if (funseq_opt) {
					signal_sum = this_var_signal;
				}
				break;
			} else {
				if (vlength > 1) { // vlength does not fall below 1
					vlength = vlength/2;
				}
				if (!(gpuCmpIntervals(this_var_chr, this_var_end, this_ann_chr, this_ann_start, this_ann_end))) {
					if (vthis == 0) { // Don't read off the end of the array
						// int_variants = 0;
						if (funseq_opt) {
							signal_sum = 0.0;
						}
						break;
					} else {
						// Take the smaller half
						vthis = vthis-vlength;
						prev_is_greater = is_greater;
						is_greater = 0;
					}
				} else {
					if (vthis == vlength_const - 1) { // Don't read off the end of the array
						// int_variants = 0;
						if (funseq_opt) {
							signal_sum = 0.0;
						}
						break;
					} else {
						// Take the larger half
						vthis = vthis+vlength;
						prev_is_greater = is_greater;
						is_greater = 1;
					}
				}
				if (vlength == 1 && ((prev_is_greater == 1 && is_greater == 0) || (prev_is_greater == 0 && is_greater == 1))) { // No intersection
					// int_variants = 0;
					if (funseq_opt) {
						signal_sum = 0.0;
					}
					break;
				}
			}
		}
		
		if (left_pointer == -1) { // No intersection
			continue;
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
			// this_var_start = gpu_var_start[vthis];
			this_var_end = gpu_var_end[vthis];
			
			if (funseq_opt) {
				this_var_signal = gpu_var_signal[vthis];
			}
		
			// Search for intersecting variants bidirectionally
			while (this_var_chr == this_ann_chr && this_ann_start < this_var_end && this_var_end <= this_ann_end) {
				// int_variants++;
				left_pointer--;
				
				if (funseq_opt) {
					signal_sum += this_var_signal;
				}
				
				if (left_pointer == 0) { // Don't read off the leading edge
					break;
				}
				
				vthis--;
				this_var_chr = gpu_var_chr[vthis];
				// this_var_start = gpu_var_start[vthis];
				this_var_end = gpu_var_end[vthis];
				
				if (funseq_opt) {
					this_var_signal = gpu_var_signal[vthis];
				}
			}
		}
		
		if (v_anchor != vlength_const-1) {
			vthis = v_anchor;
			vthis++;
		
			// Unpack current variant
			this_var_chr = gpu_var_chr[vthis];
			// this_var_start = gpu_var_start[vthis];
			this_var_end = gpu_var_end[vthis];
			
			if (funseq_opt) {
				this_var_signal = gpu_var_signal[vthis];
			}
		
			// Search for intersecting variants bidirectionally
			while (this_var_chr == this_ann_chr && this_ann_start < this_var_end && this_var_end <= this_ann_end) {
				// int_variants++;
				right_pointer++;
				
				if (funseq_opt) {
					signal_sum += this_var_signal;
				}
				
				if (right_pointer == vlength_const-1) { // Don't read off the trailing edge
					break;
				}
				
				vthis++;
				this_var_chr = gpu_var_chr[vthis];
				// this_var_start = gpu_var_start[vthis];
				this_var_end = gpu_var_end[vthis];
				
				if (funseq_opt) {
					this_var_signal = gpu_var_signal[vthis];
				}
			}
		}
		
		// DEBUG
		// printf("Find random bins\n");
		// printf("<-- %d, %d, %d -->\n", this_ann_chr, this_ann_start, this_ann_end);
		// printf("int_variants: %d\n", int_variants);
		
		// Trimer option: do we need the reference genome context?
		int trimer = (*gpu_trimer);
		
		// curandState *d_state;
		// d_state = (curandState *)malloc(sizeof(curandState));
		int tid = threadIdx.x + blockIdx.x * blockDim.x;
		curand_init(65536, tid, 0, d_state[tid]);
		
		// Variant processing loop
		for (unsigned int k = left_pointer; k <= right_pointer; k++) {
			int new_index;
			// int pos_array[window_width];
			int pos_array_pointer = 0;
			
			// Which array do we look in?
			char *this_gpu_fasta;
			if (this_var_chr == 1) {
				this_gpu_fasta = gpu_fasta_1;
			} else if (this_var_chr == 2) {
				this_gpu_fasta = gpu_fasta_2;
			} else if (this_var_chr == 3) {
				this_gpu_fasta = gpu_fasta_3;
			} else if (this_var_chr == 4) {
				this_gpu_fasta = gpu_fasta_4;
			} else if (this_var_chr == 5) {
				this_gpu_fasta = gpu_fasta_5;
			} else if (this_var_chr == 6) {
				this_gpu_fasta = gpu_fasta_6;
			} else if (this_var_chr == 7) {
				this_gpu_fasta = gpu_fasta_7;
			} else if (this_var_chr == 8) {
				this_gpu_fasta = gpu_fasta_8;
			} else if (this_var_chr == 9) {
				this_gpu_fasta = gpu_fasta_9;
			} else if (this_var_chr == 10) {
				this_gpu_fasta = gpu_fasta_10;
			} else if (this_var_chr == 11) {
				this_gpu_fasta = gpu_fasta_11;
			} else if (this_var_chr == 12) {
				this_gpu_fasta = gpu_fasta_12;
			} else if (this_var_chr == 13) {
				this_gpu_fasta = gpu_fasta_13;
			} else if (this_var_chr == 14) {
				this_gpu_fasta = gpu_fasta_14;
			} else if (this_var_chr == 15) {
				this_gpu_fasta = gpu_fasta_15;
			} else if (this_var_chr == 16) {
				this_gpu_fasta = gpu_fasta_16;
			} else if (this_var_chr == 17) {
				this_gpu_fasta = gpu_fasta_17;
			} else if (this_var_chr == 18) {
				this_gpu_fasta = gpu_fasta_18;
			} else if (this_var_chr == 19) {
				this_gpu_fasta = gpu_fasta_19;
			} else if (this_var_chr == 20) {
				this_gpu_fasta = gpu_fasta_20;
			} else if (this_var_chr == 21) {
				this_gpu_fasta = gpu_fasta_21;
			} else if (this_var_chr == 22) {
				this_gpu_fasta = gpu_fasta_22;
			} else if (this_var_chr == 23) {
				this_gpu_fasta = gpu_fasta_23;
			} else if (this_var_chr == 24) {
				this_gpu_fasta = gpu_fasta_24;
			} else if (this_var_chr == 25) {
				this_gpu_fasta = gpu_fasta_25;
			}
			
			if (trimer) {
			
				int this_var_chr = gpu_var_chr[k];
				int this_var_end = gpu_var_end[k];
				
				char cur_nt1 = toupper(this_gpu_fasta[this_var_end-2]); // 0-based index
				char cur_nt2 = toupper(this_gpu_fasta[this_var_end-1]); // 0-based index
				char cur_nt3 = toupper(this_gpu_fasta[this_var_end]); // 0-based index
				
				// If there is an N in this string, we skip this variant
				if (cur_nt1 == 'N' || cur_nt2 == 'N' || cur_nt3 == 'N') {
					gpu_permuted_chr[k] = 26; // Sentinel value
					gpu_permuted_end[k] = 0;
					continue;
				}
				
				// Build the index list
				for (int l = this_ann_start+1; l <= this_ann_end; l++) { // 1-based index
					
					// Don't read in characters if it will read off either end
					if (l == 1 || l == hg19_coor[this_var_chr-1]) {
						continue;
					}
					
					char nt1 = toupper(this_gpu_fasta[l-2]); // 0-based index
					char nt2 = toupper(this_gpu_fasta[l-1]); // 0-based index
					char nt3 = toupper(this_gpu_fasta[l]); // 0-based index
					
					if (nt1 == cur_nt1 && nt2 == cur_nt2 && nt3 == cur_nt3 && l != this_var_end) {
						pos_array[pos_array_pointer] = l; // 1-based index
						pos_array_pointer++;
					}
				}
				
				// If no positions are available, skip this variant
				if (pos_array_pointer == 0) {
					gpu_permuted_chr[k] = 26; // Sentinel value
					gpu_permuted_end[k] = 0;
					continue;
				}
				
				// Pick new position
				float this_rand = curand_uniform(d_state[tid]);
				new_index = (int)ceil(this_rand*pos_array_pointer)-1; // Selection in interval [0,pos_array_length-1]
			} else {
				do {
					float this_rand = curand_uniform(d_state[tid]);
					new_index = this_ann_start + (int)ceil(this_rand*(this_ann_end-this_ann_start)); // Selection in interval [1,(this_ann_end-this_ann_start)], added to this_ann_start
				} while (this_gpu_fasta[new_index-1] == 'N');
			}
			
			// Send the results to the output array
			gpu_permuted_chr[k] = this_var_chr;
			
			if (trimer) {
				gpu_permuted_end[k] = pos_array[new_index];
			} else {
				gpu_permuted_end[k] = new_index;
			}
		}
	}
}

__global__ void apportionWork(int* gpu_var_chr, int* gpu_var_end, 
	double* gpu_var_signal, int* gpu_ann_chr, int* gpu_ann_start, int* gpu_ann_end, 
	int* gpu_var_arr_length, int* gpu_ann_arr_length, int* gpu_trimer, 
	int *gpu_permuted_chr, int *gpu_permuted_end, 
	int *gpu_wg_switch, curandState **d_state, char *gpu_fasta_1, char *gpu_fasta_2, 
	char *gpu_fasta_3, char *gpu_fasta_4, char *gpu_fasta_5, char *gpu_fasta_6, 
	char *gpu_fasta_7, char *gpu_fasta_8, char *gpu_fasta_9, char *gpu_fasta_10, 
	char *gpu_fasta_11, char *gpu_fasta_12, char *gpu_fasta_13, char *gpu_fasta_14, 
	char *gpu_fasta_15, char *gpu_fasta_16, char *gpu_fasta_17, char *gpu_fasta_18, 
	char *gpu_fasta_19, char *gpu_fasta_20, char *gpu_fasta_21, char *gpu_fasta_22, 
	char *gpu_fasta_23, char *gpu_fasta_24, char *gpu_fasta_25, int *pos_array) {
// __global__ void apportionWork() {

	// DEBUG
	// printf("Running thread\n");
	// *test_int_gpu = 247;
	
	// Which thread am I?
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	int total_threads = 24*128;
	
	// DEBUG
// 	if (tid == 0) {
// 		printf("%f\n", gpu_signal_pvalues[13894]);
// 	}
// 	printf("tid %d: %d\n", tid, *gpu_ann_arr_length);
// 	printf("tid %d: %d\n", tid, total_threads);
	
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
		// printf("Thread ID: %d; start index: %d; end index: %d\n", tid, start, end);
// 		if (tid > 450) {
// 			return;
// 		}
		
		intersection_kernel(start, end, gpu_var_chr, gpu_var_end, gpu_var_signal, 
			gpu_ann_chr, gpu_ann_start, gpu_ann_end, gpu_var_arr_length, gpu_trimer, 
			gpu_permuted_chr, gpu_permuted_end, gpu_wg_switch, 
			d_state, gpu_fasta_1, 
			gpu_fasta_2, gpu_fasta_3, gpu_fasta_4, 
			gpu_fasta_5, gpu_fasta_6, gpu_fasta_7, gpu_fasta_8, gpu_fasta_9, gpu_fasta_10, 
			gpu_fasta_11, gpu_fasta_12, gpu_fasta_13, gpu_fasta_14, gpu_fasta_15, gpu_fasta_16, 
			gpu_fasta_17, gpu_fasta_18, gpu_fasta_19, gpu_fasta_20, gpu_fasta_21, gpu_fasta_22, 
			gpu_fasta_23, gpu_fasta_24, gpu_fasta_25, pos_array);
	} else {
		start = tid;
		end = tid;
		
		// DEBUG: Print the thread ID, start index, and end index
		// printf("Thread ID: %d; start index: %d; end index: %d\n", tid, start, end);
		
		if (tid < length) {
			intersection_kernel(start, end, gpu_var_chr, gpu_var_end, gpu_var_signal, 
				gpu_ann_chr, gpu_ann_start, gpu_ann_end, gpu_var_arr_length, gpu_trimer, 
				gpu_permuted_chr, gpu_permuted_end, gpu_wg_switch, 
				d_state, gpu_fasta_1, 
				gpu_fasta_2, gpu_fasta_3, gpu_fasta_4, 
				gpu_fasta_5, gpu_fasta_6, gpu_fasta_7, gpu_fasta_8, gpu_fasta_9, gpu_fasta_10, 
				gpu_fasta_11, gpu_fasta_12, gpu_fasta_13, gpu_fasta_14, gpu_fasta_15, gpu_fasta_16, 
				gpu_fasta_17, gpu_fasta_18, gpu_fasta_19, gpu_fasta_20, gpu_fasta_21, gpu_fasta_22, 
				gpu_fasta_23, gpu_fasta_24, gpu_fasta_25, pos_array);
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
	
	// Is the trinucleotide preservation option enabled?
	// 0 for no preservation
	// 3 for trinucleotide preservation
	// 5 for pentanucleotide preservation
	int trimer;
	
	// Number of permuted variant datasets to create
	int num_permutations;
	
	// Width of the window in which we permute variants
	int window_width;
	
	// Minimum width allowed of the variant permutation bins
	int min_width;
	
	// File with prohibited coordinates
	// Expected format: tab(chr, start, end, ...)
	string prohibited_file;
	
	// Directory with wg FASTA files
	string fasta_dir;
	
	// File with single nucleotide variants
	// Expected format: tab(chr, start, end, ...)
	string vfile;
	
	// Directory with the output files
	// Format: tab(chr, start, end, any extra columns...)
	string outdir;
	
	// Option to specify whether to calculate wg signal scores on the permutations
	// 'y': Compute wg signal scores alongside the MOAT results
	// 'n': Do not compute wg signal scores
	char funseq_opt;
	
	// WG signal file to use. Must be bigWig format.
	string signal_file;
	
	if (argc != 10 && argc != 11) {
		fprintf(stderr, "Usage: moat_v_gpu [3mer preservation option (y/n)] [# permuted datasets] [permutation window width] [min width] [prohibited regions file] [FASTA dir] [variant file] [output directory] [wg signal option (o/p/n)] [wg signal file (optional)]. Exiting.\n");
		return 1;
	} else {

		if (argv[1][0] == 'y') {
			trimer = 3;
		} else if (argv[1][0] == 'n') {
			trimer = 0;
		} else {
			fprintf(stderr, "Invalid option for 3mer preservation option: \'%c\'. Must be either \'y\' or \'n\'. Exiting.\n", argv[1][0]);
			return 1;
		}

		num_permutations = atoi(argv[2]);
		window_width = atoi(argv[3]);
		min_width = atoi(argv[4]);
		prohibited_file = string(argv[5]);
		fasta_dir = string(argv[6]);
		vfile = string(argv[7]);
		outdir = string(argv[8]);
		funseq_opt = argv[9][0];
		
		if (funseq_opt != 'y' && funseq_opt != 'n') {
			fprintf(stderr, "Invalid option for wg signal option: \'%c\'. Must be either \'y\' or \'n\'. Exiting.\n", funseq_opt);
			return 1;
		}
		
		if (argc == 11) {
			signal_file = string(argv[10]);
		}
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
	
	struct stat pbuf;
	if (stat(prohibited_file.c_str(), &pbuf)) { // Report the error and exit
		fprintf(stderr, "Error trying to stat %s: %s\n", prohibited_file.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
// 	if (pbuf.st_size == 0) {
// 		fprintf(stderr, "Error: Prohibited regions file cannot be empty. Exiting.\n");
// 		return 1;
// 	}
	
	// Check that the FASTA directory is a valid path
	struct stat fbuf;
	if (stat(fasta_dir.c_str(), &fbuf)) { // Report the error and exit
		fprintf(stderr, "Error trying to stat %s: %s\n", fasta_dir.c_str(), strerror(errno));
		return 1;
	}
	
	// Check that the outdir is a valid path
	struct stat obuf;
	if (stat(outdir.c_str(), &obuf)) { // Make the directory ourselves
		string fix_outdir = "mkdir -p " + outdir;
		system(fix_outdir.c_str());
	}
	
	// Check bigWigAverageOverBed and wg signal file in wg signal score mode
	if (funseq_opt != 'n') {
		// Verify that bigWigAverageOverBed is in the same directory as this program
		struct stat avgbuf;
		char avgoverbed_cstr[] = "./bigWigAverageOverBed";
		if (stat(avgoverbed_cstr, &avgbuf)) {
			fprintf(stderr, "Error: bigWigAverageOverBed is not in the same directory as \"moat\". Exiting.\n");
			return 1;
		}
		
		// Verify wg signal file
		struct stat databuf;
		if (stat(signal_file.c_str(), &databuf)) { // Report the error and exit
			fprintf(stderr, "Error trying to stat %s: %s\n", signal_file.c_str(), strerror(errno));
			return 1;
		}
		// Check that the file is not empty
		if (databuf.st_size == 0) {
			fprintf(stderr, "Error: WG signal file cannot be empty. Exiting.\n");
			return 1;
		}
	}
	
	// Incompatible options
	if (funseq_opt != 'n' && signal_file.empty()) {
		fprintf(stderr, "Error: Requested use of wg signal scores without specifying ");
		fprintf(stderr, "wg signal file. Exiting.\n");
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
	// Save all columns
	char linebuf[STRSIZE];
	FILE *vfile_ptr = fopen(vfile.c_str(), "r");
	while (fgets(linebuf, STRSIZE, vfile_ptr) != NULL) {
		string line = string(linebuf);
		
		// DEBUG
		// printf("%s\n", line.c_str());
		
		// Extract chromosome, start, end, ref, alt, and any extra columns from line
		vector<string> vec;
		size_t ws_index = 0;
		while (ws_index != string::npos) {
			ws_index = line.find_first_of("\t\n");
			string in = line.substr(0, ws_index);
			// printf("%s\n", in.c_str()); // DEBUG
			vec.push_back(in);
			if (ws_index != string::npos) {
				line = line.substr(ws_index+1);
			}
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
	sort(prohibited_regions.begin(), prohibited_regions.end(), cmpIntervals);
	
	// Merge prohibited regions
	prohibited_regions = merge_intervals(prohibited_regions);
	
	// hg19 coordinates
	map<string,int> hg19_coor;
	
	hg19_coor["chr1"] = 249250621;
	hg19_coor["chr2"] = 243199373;
	hg19_coor["chr3"] =	198022430;
	hg19_coor["chr4"] =	191154276;
	hg19_coor["chr5"] =	180915260;
	hg19_coor["chr6"] =	171115067;
	hg19_coor["chr7"] =	159138663;
	hg19_coor["chr8"] =	146364022;
	hg19_coor["chr9"] =	141213431;
	hg19_coor["chr10"] = 135534747;
	hg19_coor["chr11"] = 135006516;
	hg19_coor["chr12"] = 133851895;
	hg19_coor["chr13"] = 115169878;
	hg19_coor["chr14"] = 107349540;
	hg19_coor["chr15"] = 102531392;
	hg19_coor["chr16"] = 90354753;
	hg19_coor["chr17"] = 81195210;
	hg19_coor["chr18"] = 78077248;
	hg19_coor["chr19"] = 59128983;
	hg19_coor["chr20"] = 63025520;
	hg19_coor["chr21"] = 48129895;
	hg19_coor["chr22"] = 51304566;
	hg19_coor["chrX"] = 155270560;
	hg19_coor["chrY"] = 59373566;
	hg19_coor["chrM"] = 16571;
	
	for (int i = 1; i <= 25; i++) {
		string chr = int2chr(i);
		int left = 1;
		while ((left + window_radius) < hg19_coor[chr]) {
			vector<string> vec;
			vec.push_back(chr);
			
			char left_str[STRSIZE];
			sprintf(left_str, "%d", left);
			vec.push_back(string(left_str));
			
			char right_str[STRSIZE];
			sprintf(right_str, "%d", left + window_radius);
			vec.push_back(string(right_str));
			
			ann_array.push_back(vec);
			
			left += window_radius;
		}
		
		// Fencepost
		vector<string> vec;
		vec.push_back(chr);
		
		char left_str[STRSIZE];
		sprintf(left_str, "%d", left);
		vec.push_back(string(left_str));
		
		char right_str[STRSIZE];
		sprintf(right_str, "%d", hg19_coor[chr]);
		vec.push_back(string(right_str));
		
		ann_array.push_back(vec);
	}
	
	// Remove annotations that intersect the blacklist regions
	// vector<vector<string> > ann_array_new;
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		vector<vector<string> > int_intervals = intersecting_intervals(prohibited_regions, ann_array[i]);
		sort(int_intervals.begin(), int_intervals.end(), cmpIntervals);
		vector<vector<string> > allowed_regions = subtract_intervals(ann_array[i], int_intervals);
		
		if (allowed_regions.size() > 0) {
			// Easy to replace the first one
			ann_array[i][1] = allowed_regions[0][1];
			ann_array[i][2] = allowed_regions[0][2];
			
			for (unsigned int j = 1; j < allowed_regions.size(); j++) {
				ann_array.push_back(allowed_regions[j]);
			}
		} else { // No allowed regions, mark for removal
			ann_array[i][0] = "chrNo";
		}
	}
	
	sort(ann_array.begin(), ann_array.end(), cmpIntervals);
	
	// Remove those marked for deletion
	while (ann_array[ann_array.size()-1][0] == "chrNo") {
		ann_array.erase(ann_array.end());
	}
	
		// Detect those below the min threshold and join with a contiguous neighbor
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		// Unpack
		string cur_chr = ann_array[i][0];
		int cur_start = atoi(ann_array[i][1].c_str());
		int cur_end = atoi(ann_array[i][2].c_str());
		int width = cur_end - cur_start;
		
		if (width < min_width) {
			if (i != 0 && cur_chr == ann_array[i-1][0] && atoi(ann_array[i-1][2].c_str()) == cur_start) {
				// Merge
				ann_array[i-1][2] = ann_array[i][2];
			} else if (i != ann_array.size()-1 && cur_chr == ann_array[i+1][0] && atoi(ann_array[i+1][1].c_str()) == cur_end) {
				// Merge
				ann_array[i+1][1] = ann_array[i][1];
			}
			// In any case, remove it
			ann_array[i][0] = "chrNo";
		}
	}
	
	// Removal operations again
	sort(ann_array.begin(), ann_array.end(), cmpIntervals);
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
	
		// WG SIGNAL SCORE CODE
		
// 	vector<double> signal_scores;
// 	// vector<int> signal_overcount (ann_array.size(), 0);
// 	if (funseq_opt != 'n') {
// 		
// 		// Retrieve current working directory for temporary output
// 		string funseq_outdir = exec("pwd");
// 		funseq_outdir.erase(funseq_outdir.find_last_not_of(" \n\r\t")+1);
// 		funseq_outdir += "/tmp";
// 		
// 		// Verify that temporary output directory exists, or create it if it doesn't
// 		string command = "mkdir -p " + funseq_outdir;
// 		system(command.c_str());
// 		
// 		// Produce an input file for bigWigAverageOverBed in the temporary folder
// 		string avg_infile = funseq_outdir + "/" + "avg_infile.txt";
// 		string avg_outfile = funseq_outdir + "/" + "avg_outfile.txt";
// 		int regnum = 1;
// 		FILE *avg_infile_ptr = fopen(avg_infile.c_str(), "w");
// 		for (unsigned int i = 0; i < var_array.size(); i++) {
// 			char regnum_cstr[STRSIZE];
// 			sprintf(regnum_cstr, "%d", regnum);
// 			string regnum_str = "reg" + string(regnum_cstr);
// 			fprintf(avg_infile_ptr, "%s\t%s\t%s\t%s\n", var_array[i][0].c_str(), var_array[i][1].c_str(), var_array[i][2].c_str(), regnum_str.c_str());
// 			regnum++;
// 		}
// 		fclose(avg_infile_ptr);
// 		
// 		// The actual command
// 		// Assumes bigWigAverageOverBed is in same directory
// 		command = "./bigWigAverageOverBed " + signal_file + " " + avg_infile + " " + avg_outfile;
// 		system(command.c_str());
// 		
// 		// Next command depends on OS
// 		command = "uname";
// 		char buf[STRSIZE];
// 		string os = "";
// 		FILE* pipe = popen(command.c_str(), "r");
// 		if (!pipe) throw runtime_error("Could not determine operating system. Exiting.\n");
// 		try {
// 			while (!feof(pipe)) {
// 				if (fgets(buf, STRSIZE, pipe) != NULL) {
// 					os += buf;
// 				}
// 			}
// 		} catch (...) {
// 			pclose(pipe);
// 			throw;
// 		}
// 		pclose(pipe);
// 		
// 		// DEBUG
// 		// printf("%s\n", os.c_str());
// 		if (os == "Darwin\n") { // OS = Mac OS X
// 			command = "sed -i .bak 's/^reg//g' " + avg_outfile;
// 		} else { // Assume Linux, or rather, that this command is compatible
// 			command = "sed -i 's/^reg//g' " + avg_outfile;
// 		}
// 		system(command.c_str());
// 		
// 		string avg_outfile_sorted = funseq_outdir + "/" + "avg_outfile_sorted.txt";
// 		
// 		command = "sort -n -k 1,1 " + avg_outfile + " > " + avg_outfile_sorted;
// 		system(command.c_str());
// 		
// 		// Collect sum of signal scores per annotation
// 		// vector<vector<string> > signal_output;
// 		
// 		// DEBUG
// 		// printf("Breakpoint 3a\n");
// 		
// 		// Index to track where we are in the var_array
// 		unsigned int v_index = 0;
// 		
// 		// Read the output into memory
// 		FILE *avg_outfile_ptr = fopen(avg_outfile_sorted.c_str(), "r");
// 		char linebuf_cstr[STRSIZE];
// 		while (fgets(linebuf_cstr, STRSIZE-1, avg_outfile_ptr) != NULL) {
// 			
// 			string linebuf = string(linebuf_cstr);
// 			int col_index = 0;
// 			while (col_index < 5) {
// 				unsigned int pos = linebuf.find_first_of("\t");
// 				linebuf = linebuf.substr(pos+1);
// 				col_index++;
// 			}
// 			
// 			// Now linebuf has the value we're looking for. Put it in the funseq_scores vector.
// // 			double signal_score;
// // 			sscanf(linebuf.c_str(), "%lf", &signal_score);
// 			
// // 			vector<string> temp;
// // 			temp.push_back(var_array[v_index][0]);
// // 			temp.push_back(var_array[v_index][1]);
// // 			temp.push_back(var_array[v_index][2]);
// 			var_array[v_index].push_back(linebuf);
// 			// signal_output.push_back(temp);
// 			v_index++;
// 		}
// 		if (!(feof(avg_outfile_ptr)) && ferror(avg_outfile_ptr)) { // This is an error
// 			char preamble[STRSIZE];
// 			sprintf(preamble, "There was an error reading from %s", avg_outfile_sorted.c_str());
// 			perror(preamble);
// 			return 1;
// 		}
// 		fclose(avg_outfile_ptr);
// 		
// 		// Clean up temporary folder
// 		string rm_com = "rm -rf " + funseq_outdir;
// 		system(rm_com.c_str());
// 		// DEBUG remove these lines
// 		
// 		// Sort
// 		// sort(signal_output.begin(), signal_output.end(), cmpIntervals);
// 		
// 		// DEBUG
// 		// printf("Breakpoint 3b\n");
// 		
// 		// Gather up and sum the Funseq values over each annotation
// 		unsigned int signal_var_pointer = 0;
// 		for (unsigned int i = 0; i < ann_array.size(); i++) {
// 			pair<unsigned int,unsigned int> range = intersecting_variants(var_array, ann_array[i], signal_var_pointer);
// 			signal_var_pointer = range.first;
// 			double signal_sum = 0.0;
// 			
// 			for (unsigned int j = range.first; j < range.second; j++) {
// 				double score = atof(var_array[j][3].c_str());
// 				signal_sum += score;
// 			}
// 			signal_scores.push_back(signal_sum);
// 		}
// 	}

	// END WG SIGNAL SCORE CODE
	
	// Preimport the FASTA genome
	string chr_nt[25];
	if (trimer) {
		
		for (int i = 1; i <= 25; i++) {
			string filename = fasta_dir + "/" + int2chr(i) + ".fa";
			FILE *fasta_ptr = fopen(filename.c_str(), "r");
			
			int first = 1;
			char linebuf_cstr[STRSIZE];
			while (fgets(linebuf_cstr, STRSIZE, fasta_ptr) != NULL) {
				string linebuf = string(linebuf_cstr);
				linebuf.erase(linebuf.find_last_not_of(" \n\r\t")+1);
				if (first) {
					first = 0;
					continue;
				}
				chr_nt[i-1] += linebuf;
			}
			// Check feof of fasta_ptr
			if (feof(fasta_ptr)) { // We're good
				fclose(fasta_ptr);
			} else { // It's an error
				char errstring[STRSIZE];
				sprintf(errstring, "Error reading from %s", filename.c_str());
				perror(errstring);
				return 1;
			}
		}
	}
	
	// Can malloc free the prohibited regions
	prohibited_regions.clear();
	
	// Length of variant array
	int var_arr_length = var_array.size();
	
	// Turn the vectors into int arrays for CUDA
	// Var array
	int *var_chr = (int*)malloc(var_arr_length*sizeof(int));
	// int *var_start = (int*)malloc(var_arr_length*sizeof(int));
	int *var_end = (int*)malloc(var_arr_length*sizeof(int));
	double *var_signal;
	
	if (funseq_opt == 'p') {
		var_signal = (double*)malloc(var_arr_length*sizeof(double));
	}
	
	// Lengths of annotation array
	int ann_arr_length = ann_array.size();
	
	// Ann array
	int *ann_chr = (int*)malloc(ann_arr_length*sizeof(int));
	int *ann_start = (int*)malloc(ann_arr_length*sizeof(int));
	int *ann_end = (int*)malloc(ann_arr_length*sizeof(int));
	
	// Variant array processing
	for (unsigned int i = 0; i < var_array.size(); i++) {
		// Unpack the current variant
		string cur_var_chr = var_array[i][0];
		// string cur_var_start = var_array[i][1];
		string cur_var_end = var_array[i][2];
		string cur_var_signal;
		
		if (funseq_opt == 'p') {
			cur_var_signal = var_array[i][3];
		}
		
		int cur_var_chr_int = chr2int(cur_var_chr);
		
		// var_arr_length++;
		
		// var_chr = (int*)realloc(var_chr, var_arr_length*sizeof(int));
		var_chr[i] = cur_var_chr_int;
		
		// var_start = (int*)realloc(var_start, var_arr_length*sizeof(int));
		// var_start[i] = atoi(cur_var_start.c_str());
		
		// var_end = (int*)realloc(var_end, var_arr_length*sizeof(int));
		var_end[i] = atoi(cur_var_end.c_str());
		
		if (funseq_opt == 'p') {
			var_signal[i] = atof(cur_var_signal.c_str());
		}
	}
	
	// Can malloc free the variant array
	// var_array.clear();
	
	// Annotation array processing
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		// Unpack the current variant
		string cur_ann_chr = ann_array[i][0];
		string cur_ann_start = ann_array[i][1];
		string cur_ann_end = ann_array[i][2];
		
		int cur_ann_chr_int = chr2int(cur_ann_chr);
		
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
	// int *gpu_var_start;
	int *gpu_var_end;
	double *gpu_var_signal;
	
	int *gpu_ann_chr;
	int *gpu_ann_start;
	int *gpu_ann_end;
	
	int *gpu_var_arr_length;
	int *gpu_ann_arr_length;
	
	int *gpu_trimer;
	int *gpu_num_permutations;
	
	int *gpu_permuted_chr;
	int *gpu_permuted_end;
	
	// Instantiate 25 char arrays on GPU
	char *fasta_1 = chr_nt[0].c_str();
	char *fasta_2 = chr_nt[1].c_str();
	char *fasta_3 = chr_nt[2].c_str();
	char *fasta_4 = chr_nt[3].c_str();
	char *fasta_5 = chr_nt[4].c_str();
	char *fasta_6 = chr_nt[5].c_str();
	char *fasta_7 = chr_nt[6].c_str();
	char *fasta_8 = chr_nt[7].c_str();
	char *fasta_9 = chr_nt[8].c_str();
	char *fasta_10 = chr_nt[9].c_str();
	char *fasta_11 = chr_nt[10].c_str();
	char *fasta_12 = chr_nt[11].c_str();
	char *fasta_13 = chr_nt[12].c_str();
	char *fasta_14 = chr_nt[13].c_str();
	char *fasta_15 = chr_nt[14].c_str();
	char *fasta_16 = chr_nt[15].c_str();
	char *fasta_17 = chr_nt[16].c_str();
	char *fasta_18 = chr_nt[17].c_str();
	char *fasta_19 = chr_nt[18].c_str();
	char *fasta_20 = chr_nt[19].c_str();
	char *fasta_21 = chr_nt[20].c_str();
	char *fasta_22 = chr_nt[21].c_str();
	char *fasta_23 = chr_nt[22].c_str();
	char *fasta_24 = chr_nt[23].c_str();
	char *fasta_25 = chr_nt[24].c_str();
	
	char *gpu_fasta_1;
	char *gpu_fasta_2;
	char *gpu_fasta_3;
	char *gpu_fasta_4;
	char *gpu_fasta_5;
	char *gpu_fasta_6;
	char *gpu_fasta_7;
	char *gpu_fasta_8;
	char *gpu_fasta_9;
	char *gpu_fasta_10;
	char *gpu_fasta_11;
	char *gpu_fasta_12;
	char *gpu_fasta_13;
	char *gpu_fasta_14;
	char *gpu_fasta_15;
	char *gpu_fasta_16;
	char *gpu_fasta_17;
	char *gpu_fasta_18;
	char *gpu_fasta_19;
	char *gpu_fasta_20;
	char *gpu_fasta_21;
	char *gpu_fasta_22;
	char *gpu_fasta_23;
	char *gpu_fasta_24;
	char *gpu_fasta_25;
		
	// DEBUG
// 	printf("Begin CUDA code\n");
// 	int *test_int_cpu;
// 	int *test_int_gpu;
// 	test_int_cpu = (int*)malloc(sizeof(int));
// 	cudaMalloc((void**)&test_int_gpu, sizeof(int));
// 	*test_int_cpu = 246;
// 	cudaMemcpy(test_int_gpu, test_int_cpu, sizeof(int), cudaMemcpyHostToDevice);

	cudaDeviceSetLimit(cudaLimitMallocHeapSize, 8589934592);
	GPUerrchk(cudaPeekAtLastError());
	
	cudaMalloc((void**)&gpu_var_chr, var_arr_length*sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
// 	cudaMalloc((void**)&gpu_var_start, var_arr_length*sizeof(int));
// 	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_var_end, var_arr_length*sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	
	if (funseq_opt == 'p') {
		cudaMalloc((void**)&gpu_var_signal, var_arr_length*sizeof(double));
	}
	
	cudaMalloc((void**)&gpu_ann_chr, ann_arr_length*sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_ann_start, ann_arr_length*sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_ann_end, ann_arr_length*sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	
	cudaMalloc((void**)&gpu_var_arr_length, sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_ann_arr_length, sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	
	cudaMalloc((void**)&gpu_trimer, sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_num_permutations, sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	
	// Allocations for output
	cudaMalloc((void**)&gpu_permuted_chr, var_arr_length*sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_permuted_end, var_arr_length*sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	
	// Allocations for FASTA
	cudaMalloc((void**)&gpu_fasta_1, chr_nt[0].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_2, chr_nt[1].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_3, chr_nt[2].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_4, chr_nt[3].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_5, chr_nt[4].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_6, chr_nt[5].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_7, chr_nt[6].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_8, chr_nt[7].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_9, chr_nt[8].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_10, chr_nt[9].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_11, chr_nt[10].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_12, chr_nt[11].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_13, chr_nt[12].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_14, chr_nt[13].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_15, chr_nt[14].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_16, chr_nt[15].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_17, chr_nt[16].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_18, chr_nt[17].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_19, chr_nt[18].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_20, chr_nt[19].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_21, chr_nt[20].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_22, chr_nt[21].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_23, chr_nt[22].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_24, chr_nt[23].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	cudaMalloc((void**)&gpu_fasta_25, chr_nt[24].size()*sizeof(char));
	GPUerrchk(cudaPeekAtLastError());
	
// 	double *gpu_pvalues;
// 	cudaMalloc((void**)&gpu_pvalues, ann_arr_length*sizeof(double));
// 	GPUerrchk(cudaPeekAtLastError());
// 	
// 	double *gpu_signal_pvalues;
// 	if (funseq_opt == 'p') {
// 		cudaMalloc((void**)&gpu_signal_pvalues, ann_arr_length*sizeof(double));
// 	}
	
	cudaMemcpy(gpu_var_chr, var_chr, var_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
// 	cudaMemcpy(gpu_var_start, var_start, var_arr_length*sizeof(int), cudaMemcpyHostToDevice);
// 	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_var_end, var_end, var_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	if (funseq_opt == 'p') {
		cudaMemcpy(gpu_var_signal, var_signal, var_arr_length*sizeof(double), cudaMemcpyHostToDevice);
		GPUerrchk(cudaPeekAtLastError());
	}

	cudaMemcpy(gpu_ann_chr, ann_chr, ann_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_ann_start, ann_start, ann_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_ann_end, ann_end, ann_arr_length*sizeof(int), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	
	cudaMemcpy(gpu_var_arr_length, &var_arr_length, sizeof(int), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_ann_arr_length, &ann_arr_length, sizeof(int), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	
	cudaMemcpy(gpu_trimer, &trimer, sizeof(int), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_num_permutations, &num_permutations, sizeof(int), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	
	// cudaMemcpy FASTA
	cudaMemcpy(gpu_fasta_1, fasta_1, chr_nt[0].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_2, fasta_2, chr_nt[1].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_3, fasta_3, chr_nt[2].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_4, fasta_4, chr_nt[3].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_5, fasta_5, chr_nt[4].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_6, fasta_6, chr_nt[5].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_7, fasta_7, chr_nt[6].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_8, fasta_8, chr_nt[7].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_9, fasta_9, chr_nt[8].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_10, fasta_10, chr_nt[9].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_11, fasta_11, chr_nt[10].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_12, fasta_12, chr_nt[11].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_13, fasta_13, chr_nt[12].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_14, fasta_14, chr_nt[13].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_15, fasta_15, chr_nt[14].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_16, fasta_16, chr_nt[15].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_17, fasta_17, chr_nt[16].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_18, fasta_18, chr_nt[17].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_19, fasta_19, chr_nt[18].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_20, fasta_20, chr_nt[19].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_21, fasta_21, chr_nt[20].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_22, fasta_22, chr_nt[21].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_23, fasta_23, chr_nt[22].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_24, fasta_24, chr_nt[23].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_fasta_25, fasta_25, chr_nt[24].size()*sizeof(char), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	
	// Switch to indicate whether to do wg signal scores on the permuted data
	int wg_switch;
	if (funseq_opt == 'p') {
		wg_switch = 1;
	} else {
		wg_switch = 0;
	}
	int *gpu_wg_switch;
	cudaMalloc((void**)&gpu_wg_switch, sizeof(int));
	GPUerrchk(cudaPeekAtLastError());
	cudaMemcpy(gpu_wg_switch, &wg_switch, sizeof(int), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	
	// Set up thread blocks
	int num_blocks = 24;
	int threads_per_block = 128;
	int num_threads = num_blocks * threads_per_block;
	
	// DEBUG
	// printf("Breakpoint 4a\n");
	
	// Malloc additional variables to improve performance
	curandState **d_state;
	cudaMalloc((void**)&d_state, num_threads*sizeof(curandState*));
	GPUerrchk(cudaPeekAtLastError());
	
	curandState **d_state_b;
	d_state_b = (curandState**)malloc(num_threads*sizeof(curandState*));
// 	cudaMemcpy(d_state_b, d_state, num_threads*sizeof(curandState*), cudaMemcpyDeviceToHost);
// 	GPUerrchk(cudaPeekAtLastError()); 
	
	// DEBUG
	// printf("Breakpoint 4a-1\n");
	
	for (int i = 0; i < num_threads; i++) {
		// printf("Loop iter: %d\n", i); // DEBUG
		cudaMalloc((void**)&d_state_b[i], sizeof(curandState));
		GPUerrchk(cudaPeekAtLastError());
	}
	
	// DEBUG
	// printf("Breakpoint 4a-2\n");
	
	cudaMemcpy(d_state, d_state_b, num_threads*sizeof(curandState*), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	
	int **pos_array, **pos_array_b;
	
	cudaMalloc((void**)&pos_array, num_threads*sizeof(int *));
	GPUerrchk(cudaPeekAtLastError());
	
	pos_array_b = (int**)malloc(num_threads*sizeof(int *));
	
	for (int i = 0; i < num_threads; i++) {
		cudaMalloc((void**)&pos_array_b[i], 2*window_width*sizeof(int));
		GPUerrchk(cudaPeekAtLastError());
	}
	
	cudaMemcpy(pos_array, pos_array_b, num_threads*sizeof(int*), cudaMemcpyHostToDevice);
	GPUerrchk(cudaPeekAtLastError());
	
	// DEBUG
	// printf("Breakpoint 4b\n");
	
// 	int **upstream_start, **downstream_start, **upstream_start_b, **downstream_start_b;
// 	
// 	cudaMalloc((void**)&upstream_start, num_threads*sizeof(int *));
// 	GPUerrchk(cudaPeekAtLastError());
// 	cudaMalloc((void**)&downstream_start, num_threads*sizeof(int *));
// 	GPUerrchk(cudaPeekAtLastError());
// 	
// 	upstream_start_b = (int**)malloc(num_threads*sizeof(int *));
// 	downstream_start_b = (int**)malloc(num_threads*sizeof(int *));
// 	
// 	for (int i = 0; i < num_threads; i++) {
// 		cudaMalloc((void**)&upstream_start_b[i], (n/2)*sizeof(int));
// 		GPUerrchk(cudaPeekAtLastError());
// 		cudaMalloc((void**)&downstream_start_b[i], (n/2)*sizeof(int));
// 		GPUerrchk(cudaPeekAtLastError());
// 	}
// 	
// 	cudaMemcpy(upstream_start, upstream_start_b, num_threads*sizeof(int*), cudaMemcpyHostToDevice);
// 	GPUerrchk(cudaPeekAtLastError());
// 	cudaMemcpy(downstream_start, downstream_start_b, num_threads*sizeof(int*), cudaMemcpyHostToDevice);
// 	GPUerrchk(cudaPeekAtLastError());
// 	
// 	// DEBUG
// 	// printf("Breakpoint 4c\n");
// 	
// 	int **mergesort_array, **mergesort_array_b;
// 	
// 	cudaMalloc((void**)&mergesort_array, num_threads*sizeof(int *));
// 	GPUerrchk(cudaPeekAtLastError());
// 	
// 	mergesort_array_b = (int**)malloc(num_threads*sizeof(int *));
// 	
// 	for (int i = 0; i < num_threads; i++) {
// 		cudaMalloc((void**)&mergesort_array_b[i], (n/2)*sizeof(int));
// 		GPUerrchk(cudaPeekAtLastError());
// 	}
// 	
// 	cudaMemcpy(mergesort_array, mergesort_array_b, num_threads*sizeof(int*), cudaMemcpyHostToDevice);
// 	GPUerrchk(cudaPeekAtLastError());
	
	// DEBUG
	// var_signal, gpu_pvalues, gpu_signal_pvalues
// 	printf("This is debug print 1\n");
// 	for (unsigned int k = 0; k < var_arr_length; k++) {
// 		printf("%d: %f\n", k, var_signal[k]);
// 	}
// 	for (unsigned int k = 0; k < ann_arr_length; k++) {
// 		printf("%d 1: %f\n", k, gpu_pvalues[k]);
// 		printf("%d 2: %f\n", k, gpu_signal_pvalues[k]);
// 	}
	
	// DEBUG
	// printf("Breakpoint 5\n");
	
	// Adjust the heap size based on the size of the dataset
// 	if (ann_array.size() > 3000) {
// 		int fold = ann_array.size()/3000;
// 		int new_heap_size = fold*8000000;
// 		int new_heap_size = 800000000;
//  		cudaDeviceSetLimit(cudaLimitMallocHeapSize, new_heap_size);
// 	}
	vector<vector<string> > permuted_set;
	for (int i = 0; i < num_permutations; i++) {
		apportionWork<<<num_blocks, threads_per_block>>>(gpu_var_chr, gpu_var_end, 
			gpu_var_signal, gpu_ann_chr, gpu_ann_start, gpu_ann_end, gpu_var_arr_length, 
			gpu_ann_arr_length, gpu_trimer, gpu_permuted_chr, 
			gpu_permuted_end, gpu_wg_switch, d_state, gpu_fasta_1, gpu_fasta_2, gpu_fasta_3, gpu_fasta_4, 
			gpu_fasta_5, gpu_fasta_6, gpu_fasta_7, gpu_fasta_8, gpu_fasta_9, gpu_fasta_10, 
			gpu_fasta_11, gpu_fasta_12, gpu_fasta_13, gpu_fasta_14, gpu_fasta_15, gpu_fasta_16, 
			gpu_fasta_17, gpu_fasta_18, gpu_fasta_19, gpu_fasta_20, gpu_fasta_21, gpu_fasta_22, 
			gpu_fasta_23, gpu_fasta_24, gpu_fasta_25, pos_array);
		GPUerrchk(cudaPeekAtLastError());
	
		// Collect the output values, will end with same size as var_array
		int *permuted_chr = (int *)malloc(var_array.size()*sizeof(int));
		int *permuted_end = (int *)malloc(var_array.size()*sizeof(int));
	
		cudaMemcpy(permuted_chr, gpu_permuted_chr, var_arr_length*sizeof(int), cudaMemcpyDeviceToHost);
		GPUerrchk(cudaPeekAtLastError());
		cudaMemcpy(permuted_end, gpu_permuted_end, var_arr_length*sizeof(int), cudaMemcpyDeviceToHost);
		GPUerrchk(cudaPeekAtLastError());
		
		// Post-processing code
		for (int j = 0; j < var_arr_length; j++) {
			if (permuted_chr[j] == 26) {
				continue;
			}
			string cur_chr = int2chr(permuted_chr[j]);
			
			char cur_start_cstr[STRSIZE];
			sprintf(cur_start_cstr, "%d", permuted_end[j]-1);
			string cur_start = string(cur_start_cstr);
			
			char cur_end_cstr[STRSIZE];
			sprintf(cur_end_cstr, "%d", permuted_end[j]);
			string cur_end = string(cur_end_cstr);
			
			vector<string> vec;
			vec.push_back(cur_chr);
			vec.push_back(cur_start);
			vec.push_back(cur_end);
			
			for (unsigned int l = 3; l < var_array[j].size(); l++) {
				vec.push_back(var_array[j][l]);
			}
			
			permuted_set.push_back(vec);
		}
		
		// Output generation
		char perm_num[STRSIZE];
		sprintf(perm_num, "%d", i+1);
		
		string outfile = outdir + "/permutation_" + string(perm_num) + ".txt";
		FILE *outfile_ptr = fopen(outfile.c_str(), "w");
		for (unsigned int k = 0; k < permuted_set.size(); k++) {
			string outline = "";
			for (unsigned int l = 0; l < permuted_set[k].size(); l++) {
				outline += permuted_set[k][l];
				if (l < permuted_set[k].size()-1) {
					outline += "\t";
				} else {
					outline += "\n";
				}
			}
			
			fprintf(outfile_ptr, "%s", outline.c_str());
			// fprintf(outfile_ptr, "%s\t%s\t%s\n", permuted_set[k][0].c_str(), permuted_set[k][1].c_str(), permuted_set[k][2].c_str());
		}
		fclose(outfile_ptr);
	}
	
// 	double *signal_pvalues;
// 	
// 	if (funseq_opt == 'p') {
// 		signal_pvalues = (double *)malloc(ann_arr_length*sizeof(double));
// // 		if (signal_pvalues == NULL) {
// // 			printf("Houston, we have a problem\n");
// // 		}
// 		cudaMemcpy(signal_pvalues, gpu_signal_pvalues, ann_arr_length*sizeof(double), cudaMemcpyDeviceToHost);
// 		GPUerrchk(cudaPeekAtLastError());
// 	}
// 	// GPUerrchk(cudaPeekAtLastError());
// //  	if (gpu_pvalues == NULL) {
// //  		printf("Malloc fail!\n");
// //  	}
// // 	return 0;
// // 	
// // 	if (ann_arr_length < block) {
// // 		cudaMemcpy(pvalues, gpu_pvalues, ann_arr_length*sizeof(double), cudaMemcpyDeviceToHost);
// // 		GPUerrchk(cudaPeekAtLastError());
// // 	} else {
// // 		double *pvalues_ptr = pvalues;
// // 		double *gpu_pvalues_ptr = gpu_pvalues;
// // 		for (int k = 0; k < ann_arr_length; k += block) {
// // 		
// // 			// DEBUG
// // 			printf("k: %d\n", k);
// // 		
// // 			int copyblock;
// // 			if (k < ann_arr_length-block) {
// // 				copyblock = ann_arr_length-k;
// // 			} else {
// // 				copyblock = block;
// // 			}
// // 			cudaMemcpy(pvalues_ptr, gpu_pvalues_ptr, copyblock*sizeof(double), cudaMemcpyDeviceToHost);
// // 			GPUerrchk(cudaPeekAtLastError());
// // 			pvalues_ptr += block*sizeof(double);
// // 			gpu_pvalues_ptr += block*sizeof(double);
// // 		}
// // 	}
// 	// GPUerrchk(cudaDeviceSynchronize());
// 	
// 	// DEBUG
// 	// printf("Passed CUDA memcpy\n");
// // 	for (int i = 0; i < ann_arr_length; i++) {
// // 		printf("Pvalue %d: %f\n", i, gpu_pvalues[i]);
// // 	}
// 
// 	// DEBUG
// 	// printf("Breakpoint 7\n");
// 
// 	// Output generation
// 	FILE *outfile_ptr = fopen(outfile.c_str(), "w");
// 	for (unsigned int i = 0; i < ann_array.size(); i++) {
// 		
// 		// Unpack the annotation
// 		string cur_ann_chr = ann_array[i][0];
// 		string cur_ann_start = ann_array[i][1];
// 		string cur_ann_end = ann_array[i][2];
// 		string cur_ann_name = ann_array[i][3];
// 		
// 		// Print the output line
// 		if (funseq_opt == 'p') {
// 			fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%f\t%e\t%f\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str(), pvalues[i], signal_scores[i], signal_pvalues[i]);
// 		} else if (funseq_opt == 'o') {
// 			fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%f\t%e\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str(), pvalues[i], signal_scores[i]);
// 		} else {
// 			fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%f\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str(), pvalues[i]);
// 		}
// 	}
// 	fclose(outfile_ptr);
	
 	return 0;
}
