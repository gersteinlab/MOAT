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
#include <sstream>
#include <limits.h>
#include <stdexcept>
#include "variant_permutation_v3.h"
#include <mpi.h>

using namespace std;

#define STRSIZE 512
#define BIGSTRSIZE 10240

// The notable difference with the v2 is that this code places random bins anywhere
// in the genome
// In v3, rather than permuting the window locations, we instead permute the variant
// positions

// In v4, the variants are treated as positions on a number line that is shifted
// everywhere in the genome, with wraparound between chromosomes

// In v5, only the variants within some number of basepairs upstream and downstream
// of each annotation are shuffled with random uniform distribution

// In v6, we generate whole genome bins, subtract the prohibited regions, then 
// shuffle the variants within each bin, writing an output file for each
// permutation

// In v7, we generate whole genome bins, subtract the prohibited regions, then
// shuffle the variants within each bin. Whereas v6 chose random locations uniformly
// within each bin, v7 picks variant locations that preserve the reference nucleotide
// context

// In v8, the dinucleotide context of the reference influences the available
// locations

// In v9, the trinucleotide context of the reference influences the available
// locations

// vector<string> coor_search (long raw_rand_start, int ann_length) {
// 	
// 	// hg19 chromosome lengths
// 	int chr_lengths[24] = {248956422, 242193529, 198295559, 190214555, 181538259,
// 													170805979, 159345973, 145138636, 138394717, 133797422,
// 													135086622, 133275309, 114364328, 107043718, 101991189,
// 													90338345, 83257441, 80373285, 58617616, 64444167, 46709983,
// 													50818468, 156040895, 57227415};
// 	int i = 0;
// 	for (; i < 24; i++) {
// 		if (raw_rand_start >= (chr_lengths[i] - ann_length)) {
// 			raw_rand_start -= chr_lengths[i] - ann_length;
// 		} else {
// 			break;
// 		}
// 	}
// 	
// 	string rand_chr;
// 	if (i == 22) {
// 		rand_chr = "chrX";
// 	} else if (i == 23) {
// 		rand_chr = "chrY";
// 	} else {
// 		i++;
// 		char rand_chr_cstr[STRSIZE];
// 		sprintf(rand_chr_cstr, "%d", i);
// 		rand_chr = "chr" + string(rand_chr_cstr);
// 	}
// 	
// 	int rand_start = raw_rand_start;
// 	char rand_start_cstr[STRSIZE];
// 	sprintf(rand_start_cstr, "%d", rand_start);
// 	
// 	int rand_end = rand_start + ann_length;
// 	char rand_end_cstr[STRSIZE];
// 	sprintf(rand_end_cstr, "%d", rand_end);
// 	
// 	vector<string> vec;
// 	vec.push_back(rand_chr);
// 	vec.push_back(string(rand_start_cstr));
// 	vec.push_back(string(rand_end_cstr));
// 	return vec;
// }

// Circular list increment for the next_child pointer
int incrementNextChild(int next_child, int mpi_size) {
	next_child++;
	if (next_child == mpi_size) {
		next_child = 1;
	}
	return next_child;
}

// Refactorization of the code that turns a chromosome string into an integer
int chr2int (string chr_str) {
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
string int2chr (int chr_num) {
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

/*
 * Subroutine that merges the input intervals (3col)
 * Assumes input interval array is sorted
 */
vector<vector<string> > merge_intervals (vector<vector<string> > starting_intervals) {
	
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
 * Subroutine for calculating the mutation counts given a sorted vector of 
 * variants, and a sorted vector of annotations. Returns a vector of the same
 * length as the annotation vector, with the mutation counts for each annotation
 */
vector<int> mutation_counts (vector<vector<string> > var_array, vector<vector<string> > ann_array) {
	
	// Variables for main loop
	unsigned int var_pointer = 0;
	
	// Output vector
	vector<int> counts;
	
	// Main loop: Iterate through the annotations
	for (unsigned int i = 0; i < ann_array.size(); i++) {
	
		// Unpack the current annotation
		string cur_ann_chr = ann_array[i][0];
		string cur_ann_start = ann_array[i][1];
		string cur_ann_end = ann_array[i][2];
		string cur_ann_name = ann_array[i][3];
		
		int cur_ann_chr_num = chr2int(cur_ann_chr);
		int cur_ann_start_num = atoi(cur_ann_start.c_str());
		int cur_ann_end_num = atoi(cur_ann_end.c_str());
		
		// Find the intersecting variants for this annotation
		int target_variants = 0;
		
		// Start searching from var_pointer
		// Instantiate cur_var variables
		string cur_var_chr = var_array[var_pointer][0];
		string cur_var_start = var_array[var_pointer][1];
		string cur_var_end = var_array[var_pointer][2];
		
		int cur_var_chr_num = chr2int(cur_var_chr);
		int cur_var_start_num = atoi(cur_var_start.c_str());
		int cur_var_end_num = atoi(cur_var_end.c_str());
		
		// Now count the intersecting variants
		// var_pointer points to the "earliest" possible intersecting variant, and vpointer
		// points to the variants up until the last intersecting with the annotation
		unsigned int vpointer = var_pointer;
		
		// While vpointer does not go past the current annotation
		while (cur_var_chr_num < cur_ann_chr_num || (cur_var_chr_num == cur_ann_chr_num && cur_var_start_num < cur_ann_end_num)) {
			
			// If the current variant intersects the current annotation, increment target_variants
			if (cur_var_chr == cur_ann_chr && cur_ann_start_num <= cur_var_end_num && cur_var_start_num <= cur_ann_end_num) {
				target_variants++;
			} else { // Update the var_pointer
				if (vpointer != var_array.size()-1) {
					var_pointer = vpointer + 1;
				}
			}
			// Now update the cur_var
			if (vpointer == var_array.size()-1) {
				break;
			}
			vpointer++;
			
			cur_var_chr = var_array[vpointer][0];
			cur_var_start = var_array[vpointer][1];
			cur_var_end = var_array[vpointer][2];
			
			// Cast chr, start and end into int
			cur_var_chr_num = chr2int(cur_var_chr);
			cur_var_start_num = atoi(cur_var_start.c_str());
			cur_var_end_num = atoi(cur_var_end.c_str());
		}
		// Collect output
		counts.push_back(target_variants);
	}
	return counts;
}

// Get the next FASTA file in the sequence
// void newFastaFilehandle (FILE *fasta_ptr, string chr) {
// 	if (fasta_ptr != NULL) {
// 		fclose(fasta_ptr);
// 	}
// 	string filename = chr + ".fa";
// 	fasta_ptr = fopen(filename.c_str(), "r");
// 	
// 	// Get past the first comment line
// 	char linebuf[STRSIZE];
// 	fgets(linebuf, STRSIZE, fasta_ptr);
// }

// Retrieve a string of characters corresponding to reference basepairs in the
// target region
// string fillBuffer(fasta_ptr, &char_pointer, region_start, region_end, ) {
// 	string buffer = "";
// 	string strbuf;
// 	char linebuf[STRSIZE];
// 	while((*char_pointer) < region_end && fgets(linebuf, STRSIZE, fasta_ptr) != NULL) {
// 		int length = strlen(linebuf);
// 		strbuf = string(linebuf);
// 		if ((*char_pointer) <= region_start && (*char_pointer)+length > region_start) {
// 			buffer += strbuf.substr(region_start-(*char_pointer));
// 		} else if ((*char_pointer) > region_start && (*char_pointer)+length <= region_end) {
// 			buffer += strbuf;
// 		} else if ((*char_pointer) <= region_end && (*char_pointer)+length > region_end) {
// 			buffer += strbuf.substr(0, region_end-(*char_pointer)+1);
// 		}
// 		(*char_pointer) += length;
// 	}
// 	return buffer;
// }

/* 
 * This code takes as input a variant track and an annotation track. For each
 * element, count intersecting variants (n_r). Then, permute all variant positions 
 * and run the intersection and mutation count for the permuted dataset (n_pi).
 * Do this for a number of times, and calculate a p-value for mutation
 * burden based on how many of the permuted datasets have (n_pi >= n_r)
 */
 
int main (int argc, char* argv[]) {

	/* MPI setup portion */
	int mpi_size;
	int mpi_rank;
	MPI_Status status;
	
	// Spawn parallel processes
	MPI_Init(&argc, &argv);
	
	// Initialize MPI variables
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	
	if (mpi_rank == 0) { // Parent process: farm work to the child processes
	
		/* User-supplied arguments */
		
		// Is the 3mer/5mer preservation option enabled?
		int trimer;
	
		// Number of permuted variant datasets to create
		int num_permutations;
	
		// Radius of the window in which we permute variants
		int window_radius;
	
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
		// Format: tab(chr, start, end)
		string outdir;
		
		// Option to specify whether to calculate wg signal scores on the permutations
		// 'y': Compute wg signal scores alongside the MOAT results
		// 'n': Do not compute wg signal scores
		char funseq_opt;
		
		// WG signal file to use. Must be bigWig format.
		string signal_file;
	
		if (argc != 10 && argc != 11) {
			fprintf(stderr, "Usage: moat_v_mpi [3mer/5mer preservation option (0/3/5)] [# permuted datasets] [permutation window radius] [min width] [prohibited regions file] [FASTA dir] [variant file] [output directory] [wg signal option (y/n)] [wg signal file (optional)]. Exiting.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		} else {
		
// 			if (argv[1][0] == 'y') {
// 				trimer = true;
// 			} else if (argv[1][0] == 'n') {
// 				trimer = false;
			if (argv[1][0] != '0' && argv[1][0] != '3' && argv[1][0] != '5') {
				fprintf(stderr, "Invalid option for 3mer preservation option: \'%c\'. Must be one of 0, 3, or 5. Exiting.\n", argv[1][0]);
				MPI_Abort(MPI_COMM_WORLD, 1);
				return 1;
			}
		
			trimer = atoi(argv[1]);
			num_permutations = atoi(argv[2]);
			window_radius = atoi(argv[3]);
			min_width = atoi(argv[4]);
			prohibited_file = string(argv[5]);
			fasta_dir = string(argv[6]);
			vfile = string(argv[7]);
			outdir = string(argv[8]);
			funseq_opt = argv[9][0];
			
			if (funseq_opt != 'y' && funseq_opt != 'n') {
				fprintf(stderr, "Invalid option for wg signal option: \'%c\'. Must be either \'y\' or \'n\'. Exiting.\n", funseq_opt);
				MPI_Abort(MPI_COMM_WORLD, 1);
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
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		// Check that the file is not empty
		if (vbuf.st_size == 0) {
			fprintf(stderr, "Error: Variant file cannot be empty. Exiting.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
	
		struct stat pbuf;
		if (stat(prohibited_file.c_str(), &pbuf)) { // Report the error and exit
			fprintf(stderr, "Error trying to stat %s: %s\n", prohibited_file.c_str(), strerror(errno));
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		// Check that the file is not empty
		if (pbuf.st_size == 0) {
			fprintf(stderr, "Error: Prohibited regions file cannot be empty. Exiting.\n");
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		
		// Check that the FASTA directory is a valid path
		struct stat fbuf;
		if (stat(fasta_dir.c_str(), &fbuf)) { // Report the error and exit
			fprintf(stderr, "Error trying to stat %s: %s\n", fasta_dir.c_str(), strerror(errno));
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		
		// Check that the outdir is a valid path
		struct stat obuf;
		if (stat(outdir.c_str(), &obuf)) { // Report the error and exit
			fprintf(stderr, "Error trying to stat %s: %s\n", outdir.c_str(), strerror(errno));
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
		
		// Check bigWigAverageOverBed and Funseq data file in static mode
		if (funseq_opt != 'n') {
			// Verify that bigWigAverageOverBed is in the same directory as this program
			struct stat avgbuf;
			char avgoverbed_cstr[] = "./bigWigAverageOverBed";
			if (stat(avgoverbed_cstr, &avgbuf)) {
				fprintf(stderr, "Error: bigWigAverageOverBed is not in the same directory as \"moat\". Exiting.\n");
				MPI_Abort(MPI_COMM_WORLD, 1);
				return 1;
			}
		
			// Verify precomputed Funseq scores are in the same directory as this program
			struct stat databuf;
			if (stat(signal_file.c_str(), &databuf)) { // Report the error and exit
				fprintf(stderr, "Error trying to stat %s: %s\n", signal_file.c_str(), strerror(errno));
				MPI_Abort(MPI_COMM_WORLD, 1);
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
		// Variant array, contains variants of the format vector(chr, start, end)
		vector<vector<string> > var_array;
	
		// Annotation array, contains annotations of the format vector(chr, start, end, ann_name)
		// Will be generated from the whole genome coordinates at runtime
		vector<vector<string> > ann_array;
	
		// Prohibited regions array, contains annotations of the format vector(chr, start, end)
		vector<vector<string> > prohibited_regions;
	
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
			MPI_Abort(MPI_COMM_WORLD, 1);
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
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}
	
		// DEBUG
		printf("Breakpoint 1\n");
	
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
	
		// DEBUG - check ann_array values
// 		FILE *testfile_ptr = fopen("test-bin-code/testfile.txt", "w");
// 		for (unsigned int i = 0; i < ann_array.size(); i++) {
// 			fprintf(testfile_ptr, "%s\t%s\t%s\n", ann_array[i][0].c_str(), ann_array[i][1].c_str(), ann_array[i][2].c_str());
// 		}
// 		fclose(testfile_ptr);
// 		return 0;
	
		// Subtract the portions of each annotation that overlap prohibited regions
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
		
		// DEBUG - check ann_array values
// 		FILE *testfile_ptr = fopen("test-bin-code/testfile.txt", "w");
// 		for (unsigned int i = 0; i < ann_array.size(); i++) {
// 			fprintf(testfile_ptr, "%s\t%s\t%s\n", ann_array[i][0].c_str(), ann_array[i][1].c_str(), ann_array[i][2].c_str());
// 		}
// 		fclose(testfile_ptr);
// 		return 0;
	
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
		printf("Breakpoint 2\n");
	
		// DEBUG - check ann_array values after prohibited region subtraction
// 		FILE *testfile_ptr = fopen("test-bin-code/testfile.txt", "w");
// 		for (unsigned int i = 0; i < ann_array.size(); i++) {
// 			fprintf(testfile_ptr, "%s\t%s\t%s\n", ann_array[i][0].c_str(), ann_array[i][1].c_str(), ann_array[i][2].c_str());
// 		}
// 		fclose(testfile_ptr);
// 		return 0;
	
// 		FILE *fasta_ptr = NULL;
// 		string last_chr = "";
// 		// int char_pointer;
// 		string chr_nt;
		
		// Pointer into a circular list that indicates the next available child process
		// Has corresponding increment method that ensures the circular list property
		// int next_child = 1;
		
		// First, give the trimer boolean flag to all children
		MPI_Bcast(&trimer, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		// Second, give the FASTA directory location to all children
		int strlen = fasta_dir.size() + 1;
		char *fasta_dir_cstr = (char *)malloc((strlen)*sizeof(char));
		strcpy(fasta_dir_cstr, fasta_dir.c_str());
		MPI_Bcast(&strlen, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(fasta_dir_cstr, strlen, MPI_CHAR, 0, MPI_COMM_WORLD);
	
		/* Permute variant locations */
		srand(0);
	
		for (int i = 0; i < num_permutations; i++) {
		
			/* Signal child processes to begin */
			int permutation_flag = 1;
			MPI_Bcast(&permutation_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
			/* Begin dividing work for MPI processes */
			
			// Vector of starting indices
			vector<unsigned int> chr_starters (25,UINT_MAX);
			
			// Iterate through the chromosomes
			unsigned int chr_ann_pointer = 0;
			unsigned int chr_var_pointer = 0;
			for (int j = 1; j <= 25; j++) {
				string chr = int2chr(j);
				
				// DEBUG
				// printf("chr: %s\n", chr.c_str());
				
				// Prepare the variant data
				vector<vector<string> > chr_var_array;
				
				for (unsigned int k = chr_var_pointer; k < var_array.size(); k++) {
					if (chr2int(var_array[k][0]) == j) {
						if (chr_starters[j-1] == UINT_MAX) {
							chr_starters[j-1] = k;
						}
						chr_var_array.push_back(var_array[k]);
					} else if (chr2int(var_array[k][0]) > j) {
						break;
					} // else this_chr < j, need to catch up
					chr_var_pointer++;
				}
				
				int chr_var_array_size = chr_var_array.size();
				
				// DEBUG
				// printf("chr_var_array_size: %d\n", chr_var_array_size);
				
				// Proceed only if there are a nonzero number of variants
				if (chr_var_array_size > 0) {
				
					// Prepare the annotation data
					vector<vector<string> > chr_ann_array;
				
					for (unsigned int k = chr_ann_pointer; k < ann_array.size(); k++) {
						if (chr2int(ann_array[k][0]) == j) {
							chr_ann_array.push_back(ann_array[k]);
						} else if (chr2int(ann_array[k][0]) > j) {
							break;
						} // else this_chr < j, need to catch up
						chr_ann_pointer++;
					}
				
					int chr_ann_array_size = chr_ann_array.size();
					
					// DEBUG: check ann array sizes
					// printf("Perm: %d, chr: %s, chr_ann_array_size: %d\n", i, chr.c_str(), chr_ann_array_size);
				
					// Determine which child is available to take this job
					int available_flag;
					MPI_Recv(&available_flag, 1, MPI_INT, MPI_ANY_SOURCE, 9, MPI_COMM_WORLD, &status);
					int next_child = status.MPI_SOURCE;
				
					// Transmit the flag indicating parent sending/child receiving (tag = 5)
					// 0 = parent send, 1 = parent receive
					int flag = 0;
					MPI_Send(&flag, 1, MPI_INT, next_child, 5, MPI_COMM_WORLD);
				
					// Transmit chromosome number to the next available child
					MPI_Send(&j, 1, MPI_INT, next_child, 0, MPI_COMM_WORLD);
				
					// Package the annotation data for the child
					int *chr_ann_coor = (int *)malloc(2*chr_ann_array_size*sizeof(int));
					for (int k = 0; k < chr_ann_array_size; k++) {
						chr_ann_coor[2*k] = atoi(chr_ann_array[k][1].c_str());
						chr_ann_coor[2*k+1] = atoi(chr_ann_array[k][2].c_str());
					}
				
					// Transmit annotation data
					int chr_ann_coor_size = 2*chr_ann_array_size;
					MPI_Send(&chr_ann_coor_size, 1, MPI_INT, next_child, 1, MPI_COMM_WORLD);
					MPI_Send(chr_ann_coor, 2*chr_ann_array_size, MPI_INT, next_child, 2, MPI_COMM_WORLD);
				
					// Package the variant data for the child
					int *chr_var_coor = (int *)malloc(chr_var_array_size*sizeof(int));
					// char *chr_var_alleles = (char *)malloc(2*chr_var_array_size*(sizeof(char)));
					for (int k = 0; k < chr_var_array_size; k++) {
						chr_var_coor[k] = atoi(chr_var_array[k][2].c_str());
						// chr_var_alleles[2*k] = chr_var_array[k][3][0];
						// chr_var_alleles[2*k+1] = chr_var_array[k][4][0];
					}
				
					// Transmit variant data
					MPI_Send(&chr_var_array_size, 1, MPI_INT, next_child, 3, MPI_COMM_WORLD);
					MPI_Send(chr_var_coor, chr_var_array_size, MPI_INT, next_child, 4, MPI_COMM_WORLD);
					// MPI_Send(chr_var_alleles, 2*chr_var_array_size, MPI_CHAR, next_child, 20, MPI_COMM_WORLD);
				
					// Free dynamic memory allocations
					free(chr_ann_coor);
					free(chr_var_coor);
					// free(chr_var_alleles);
				
					// next_child = incrementNextChild(next_child, mpi_size);
				}
			}
			
			// DEBUG
			printf("Breakpoint 1\n");
			
			// NEW CODE
			// Do wg signal step on the previously created permutation
			vector<double> funseq_scores;
			if (funseq_opt != 'n') {
				if (i >= 1) {
				
					// Retrieve current working directory for temporary Funseq2 output
					string funseq_outdir = exec("pwd");
					funseq_outdir.erase(funseq_outdir.find_last_not_of(" \n\r\t")+1);
					funseq_outdir += "/tmp";
	
					// Verify that funseq output directory exists, or create it if it doesn't
					string command = "mkdir -p " + funseq_outdir;
					system(command.c_str());
					
					// Get the path to the file to use as Funseq input
					char perm_num[STRSIZE];
					sprintf(perm_num, "%d", i);
					string outfile = outdir + "/permutation_" + string(perm_num) + ".txt";
				
					// Convert outfile to absolute path, if it is not already
					char abs_path_outfile[PATH_MAX];
					errno = 0;
					realpath(outfile.c_str(), abs_path_outfile);
					if (errno) {
						fprintf(stderr, "Error resolving absolute path of permutation variant file: %s\n", strerror(errno));
						fprintf(stderr, "Exiting.\n");
						MPI_Abort(MPI_COMM_WORLD, 1);
						return 1;
					}
					
					// Produce an input file for bigWigAverageOverBed in the temporary Funseq folder
					string avg_infile = funseq_outdir + "/" + "avg_infile.txt";
					string avg_outfile = funseq_outdir + "/" + "avg_outfile.txt";
					int regnum = 1;
					FILE *avg_infile_ptr = fopen(avg_infile.c_str(), "w");
					
					// Need to read the other columns from the original outfile
					char linebuf3[STRSIZE];
					FILE *outfile_ptr = fopen(outfile.c_str(), "r");
					while (fgets(linebuf3, STRSIZE, outfile_ptr) != NULL) {
					
						string line = string(linebuf3);
						
						// Extract chromosome, start, and end (first 3 columns)
						vector<string> vec;
						for (int i = 0; i < 3; i++) {
							size_t ws_index = line.find_first_of("\t\n");
							string in = line.substr(0, ws_index);
							vec.push_back(in);
							line = line.substr(ws_index+1);
						}
						
						char regnum_cstr[STRSIZE];
						sprintf(regnum_cstr, "%d", regnum);
						string regnum_str = "reg" + string(regnum_cstr);
						
						fprintf(avg_infile_ptr, "%s\t%s\t%s\t%s\n", vec[0].c_str(), vec[1].c_str(), vec[2].c_str(), regnum_str.c_str());
						regnum++;
					}
					// Check feof of outfile
					if (feof(outfile_ptr)) { // We're good
						fclose(outfile_ptr);
					} else { // It's an error
						char errstring[STRSIZE];
						sprintf(errstring, "Error reading from %s", outfile.c_str());
						perror(errstring);
						MPI_Abort(MPI_COMM_WORLD, 1);
						return 1;
					}
					
					fclose(avg_infile_ptr);
						
					// The actual command
					// Assumes bigWigAverageOverBed is in same directory
					command = "./bigWigAverageOverBed " + signal_file + " "  + avg_infile + " " + avg_outfile;
					system(command.c_str());
			
					// Next command depends on OS
					command = "uname";
					char buf[STRSIZE];
					string os = "";
					FILE* pipe = popen(command.c_str(), "r");
					if (!pipe) throw runtime_error("Could not determine operating system. Exiting.\n");
					try {
						while (!feof(pipe)) {
							if (fgets(buf, STRSIZE, pipe) != NULL) {
								os += buf;
							}
						}
					} catch (...) {
						pclose(pipe);
						throw;
					}
					pclose(pipe);
		
					// DEBUG
					// printf("%s\n", os.c_str());
					if (os == "Darwin\n") { // OS = Mac OS X
						command = "sed -i .bak 's/^reg//g' " + avg_outfile;
					} else { // Assume Linux, or rather, that this command is compatible
						command = "sed -i 's/^reg//g' " + avg_outfile;
					}
					system(command.c_str());
					
					string avg_outfile_sorted = funseq_outdir + "/" + "avg_outfile_sorted.txt";
			
					command = "sort -n -k 1,1 " + avg_outfile + " > " + avg_outfile_sorted;
					system(command.c_str());
			
					// Read the output into memory
					FILE *avg_outfile_ptr = fopen(avg_outfile_sorted.c_str(), "r");
					char linebuf_cstr[STRSIZE];
					while (fgets(linebuf_cstr, STRSIZE-1, avg_outfile_ptr) != NULL) {
			
						string linebuf = string(linebuf_cstr);
						int col_index = 0;
						while (col_index < 5) {
							unsigned int pos = linebuf.find_first_of("\t");
							linebuf = linebuf.substr(pos+1);
							col_index++;
						}
			
						// Now linebuf has the value we're looking for. Put it in the funseq_scores vector.
						double funseq_score;
						sscanf(linebuf.c_str(), "%lf", &funseq_score);
						funseq_scores.push_back(funseq_score);
					}
					if (!(feof(avg_outfile_ptr)) && ferror(avg_outfile_ptr)) { // This is an error
						char preamble[STRSIZE];
						sprintf(preamble, "There was an error reading from %s", avg_outfile_sorted.c_str());
						perror(preamble);
						return 1;
					}
					fclose(avg_outfile_ptr);
					
					// Print the Funseq scores to a new file that will replace the old one
					string funseq_outfile = outdir + "/permutation_" + string(perm_num) + ".funseq.txt";
					FILE *funseq_outfile_ptr = fopen(funseq_outfile.c_str(), "w");
				
					// Funseq score index
					unsigned int fs_index = 0;
				
					// Need to read the other columns from the original outfile
					// linebuf3[STRSIZE];
					outfile_ptr = fopen(outfile.c_str(), "r");
					while (fgets(linebuf3, STRSIZE, outfile_ptr) != NULL) {
						string line = string(linebuf3);
					
						// DEBUG: check what's coming in from this file
						// printf("%s\n", line.c_str());
					
						// Extract chromosome, start, end, ref, and alt (first 5 columns)
						vector<string> vec;
						for (int i = 0; i < 5; i++) {
							size_t ws_index = line.find_first_of("\t\n");
							string in = line.substr(0, ws_index);
							vec.push_back(in);
							line = line.substr(ws_index+1);
						}
					
						fprintf(funseq_outfile_ptr, "%s\t%s\t%s\t%e\n", vec[0].c_str(), vec[1].c_str(), vec[2].c_str(), funseq_scores[fs_index]);
						fs_index++;
					}
					// Check feof of outfile
					if (feof(outfile_ptr)) { // We're good
						fclose(outfile_ptr);
					} else { // It's an error
						char errstring[STRSIZE];
						sprintf(errstring, "Error reading from %s", outfile.c_str());
						perror(errstring);
						MPI_Abort(MPI_COMM_WORLD, 1);
						return 1;
					}
	
	// 					for (unsigned int k = 0; k < funseq_scores.size(); k++) {
	// 						fprintf(funseq_outfile_ptr, "%s\t%s\t%s\t%s\t%s\t%e\n", permuted_set[k][0].c_str(), permuted_set[k][1].c_str(), permuted_set[k][2].c_str(), permuted_set[k][3].c_str(), permuted_set[k][4].c_str(), funseq_scores[k]);
	// 					}
					fclose(funseq_outfile_ptr);
	
					string file_switch_1 = "rm " + outfile;
					system(file_switch_1.c_str());
					string file_switch_2 = "mv " + funseq_outfile + " " + outfile;
					system(file_switch_2.c_str());
			
						// DEBUG
						// printf("Filesystem magic done\n");
			
					// Clean up temporary folder
					string rm_com = "rm -rf " + funseq_outdir;
					system(rm_com.c_str());
				}
			}
			
			// Begin receiving data and produce output
			// First flip the send/receive roles
			int counter = 0;
			while (counter < mpi_size-1) {
				int available_flag;
				MPI_Recv(&available_flag, 1, MPI_INT, MPI_ANY_SOURCE, 9, MPI_COMM_WORLD, &status);
				int next_child = status.MPI_SOURCE;
				int flag = 1;
				MPI_Send(&flag, 1, MPI_INT, next_child, 5, MPI_COMM_WORLD);
				counter++;
			}
			
			// Receiving time
			vector<vector<string> > permuted_set;
			
			// DEBUG
			printf("Breakpoint 2\n");
			
			// Count how many processes we've heard back from
			counter = 0;
			while (counter < mpi_size-1) {
				int permuted_var_coor_size;
				MPI_Recv(&permuted_var_coor_size, 1, MPI_INT, MPI_ANY_SOURCE, 6, MPI_COMM_WORLD, &status);
				int source = status.MPI_SOURCE;
				
				// DEBUG
				// printf("Permuted variant count: perm: %d, proc: %d, num: %d\n", i, source, permuted_var_coor_size);
				
				if (permuted_var_coor_size > 0) {
					int *permuted_var_coor = (int *)malloc(permuted_var_coor_size*sizeof(int));
					// char *permuted_var_alleles = (char *)malloc(permuted_var_coor_size*sizeof(char));
					MPI_Recv(permuted_var_coor, permuted_var_coor_size, MPI_INT, source, 7, MPI_COMM_WORLD, &status);
					// MPI_Recv(permuted_var_alleles, permuted_var_coor_size, MPI_CHAR, source, 21, MPI_COMM_WORLD, &status);
					
					unsigned int *permuted_var_id = (unsigned int *)malloc((permuted_var_coor_size/2)*sizeof(unsigned int));
					MPI_Recv(permuted_var_id, (permuted_var_coor_size/2), MPI_UNSIGNED, source, 23, MPI_COMM_WORLD, &status);
				
					for (int j = 0; j < permuted_var_coor_size/2; j++) {
						string this_chr = int2chr(permuted_var_coor[2*j]);
				
						char start_str[STRSIZE];
						sprintf(start_str, "%d", permuted_var_coor[2*j+1]-1);
					
						char end_str[STRSIZE];
						sprintf(end_str, "%d", permuted_var_coor[2*j+1]);
						
// 						char ref_str[STRSIZE];
// 						sprintf(ref_str, "%c", permuted_var_alleles[2*j]);
// 					
// 						char alt_str[STRSIZE];
// 						sprintf(alt_str, "%c", permuted_var_alleles[2*j+1]);
					
						vector<string> vec;
						vec.push_back(this_chr);
						vec.push_back(string(start_str));
						vec.push_back(string(end_str));
						
						for (unsigned int k = 3; k < var_array[chr_starters[chr2int(this_chr)-1]+permuted_var_id[j]].size(); k++) {
							vec.push_back(var_array[chr_starters[chr2int(this_chr)-1]+permuted_var_id[j]][k]);
						}
// 						vec.push_back(string(ref_str));
// 						vec.push_back(string(alt_str));
						permuted_set.push_back(vec);
					}

					free(permuted_var_coor);
					// free(permuted_var_alleles);
					free(permuted_var_id);
				}
				counter++;
			}
			
			// DEBUG
			printf("Breakpoint 3\n");
			
			// For DEBUG purposes, disable the sort function
			sort(permuted_set.begin(), permuted_set.end(), cmpIntervals);
			
			// DEBUG
			// printf("Breakpoint 3a\n");
		
			// Open new output file
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
			
			// DEBUG
			printf("Breakpoint 4\n");
		}
		
		/* Signal child processes to end */
		int permutation_flag = 0;
		MPI_Bcast(&permutation_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
// 		for (int j = 1; j < mpi_size; j++) {
// 			MPI_Send(&permutation_flag, 1, MPI_INT, j, 8, MPI_COMM_WORLD);
// 		}

		// Funseq fencepost
		vector<double> funseq_scores;
		if (funseq_opt != 'n') {
		
			// Retrieve current working directory for temporary Funseq2 output
			string funseq_outdir = exec("pwd");
			funseq_outdir.erase(funseq_outdir.find_last_not_of(" \n\r\t")+1);
			funseq_outdir += "/tmp";
	
			// Verify that funseq output directory exists, or create it if it doesn't
			string command = "mkdir -p " + funseq_outdir;
			system(command.c_str());
			
			// Get the path to the file to use as Funseq input
			char perm_num[STRSIZE];
			sprintf(perm_num, "%d", num_permutations);
			string outfile = outdir + "/permutation_" + string(perm_num) + ".txt";
			
			// Convert outfile to absolute path, if it is not already
			char abs_path_outfile[PATH_MAX];
			errno = 0;
			realpath(outfile.c_str(), abs_path_outfile);
			if (errno) {
				fprintf(stderr, "Error resolving absolute path of permutation variant file: %s\n", strerror(errno));
				fprintf(stderr, "Exiting.\n");
				MPI_Abort(MPI_COMM_WORLD, 1);
				return 1;
			}
			
			// Produce an input file for bigWigAverageOverBed in the temporary Funseq folder
			string avg_infile = funseq_outdir + "/" + "avg_infile.txt";
			string avg_outfile = funseq_outdir + "/" + "avg_outfile.txt";
			int regnum = 1;
			FILE *avg_infile_ptr = fopen(avg_infile.c_str(), "w");
			
			// Need to read the other columns from the original outfile
			char linebuf3[STRSIZE];
			FILE *outfile_ptr = fopen(outfile.c_str(), "r");
			while (fgets(linebuf3, STRSIZE, outfile_ptr) != NULL) {
			
				string line = string(linebuf3);
						
				// Extract chromosome, start, and end (first 3 columns)
				vector<string> vec;
				for (int i = 0; i < 3; i++) {
					size_t ws_index = line.find_first_of("\t\n");
					string in = line.substr(0, ws_index);
					vec.push_back(in);
					line = line.substr(ws_index+1);
				}
				
				char regnum_cstr[STRSIZE];
				sprintf(regnum_cstr, "%d", regnum);
				string regnum_str = "reg" + string(regnum_cstr);
				
				fprintf(avg_infile_ptr, "%s\t%s\t%s\t%s\n", vec[0].c_str(), vec[1].c_str(), vec[2].c_str(), regnum_str.c_str());
				regnum++;
			}
			// Check feof of outfile
			if (feof(outfile_ptr)) { // We're good
				fclose(outfile_ptr);
			} else { // It's an error
				char errstring[STRSIZE];
				sprintf(errstring, "Error reading from %s", outfile.c_str());
				perror(errstring);
				MPI_Abort(MPI_COMM_WORLD, 1);
				return 1;
			}
			
			fclose(avg_infile_ptr);
				
			// The actual command
			// Assumes bigWigAverageOverBed is in same directory
			command = "./bigWigAverageOverBed " + signal_file + " " + avg_infile + " " + avg_outfile;
			system(command.c_str());
	
			// Next command depends on OS
			command = "uname";
			char buf[STRSIZE];
			string os = "";
			FILE* pipe = popen(command.c_str(), "r");
			if (!pipe) throw runtime_error("Could not determine operating system. Exiting.\n");
			try {
				while (!feof(pipe)) {
					if (fgets(buf, STRSIZE, pipe) != NULL) {
						os += buf;
					}
				}
			} catch (...) {
				pclose(pipe);
				throw;
			}
			pclose(pipe);

			// DEBUG
			// printf("%s\n", os.c_str());
			if (os == "Darwin\n") { // OS = Mac OS X
				command = "sed -i .bak 's/^reg//g' " + avg_outfile;
			} else { // Assume Linux, or rather, that this command is compatible
				command = "sed -i 's/^reg//g' " + avg_outfile;
			}
			system(command.c_str());
			
			string avg_outfile_sorted = funseq_outdir + "/" + "avg_outfile_sorted.txt";
	
			command = "sort -n -k 1,1 " + avg_outfile + " > " + avg_outfile_sorted;
			system(command.c_str());
				
			// Read the output into memory
			FILE *avg_outfile_ptr = fopen(avg_outfile_sorted.c_str(), "r");
			char linebuf_cstr[STRSIZE];
			while (fgets(linebuf_cstr, STRSIZE-1, avg_outfile_ptr) != NULL) {
	
				string linebuf = string(linebuf_cstr);
				int col_index = 0;
				while (col_index < 5) {
					unsigned int pos = linebuf.find_first_of("\t");
					linebuf = linebuf.substr(pos+1);
					col_index++;
				}
	
				// Now linebuf has the value we're looking for. Put it in the funseq_scores vector.
				double funseq_score;
				sscanf(linebuf.c_str(), "%lf", &funseq_score);
				funseq_scores.push_back(funseq_score);
			}
			if (!(feof(avg_outfile_ptr)) && ferror(avg_outfile_ptr)) { // This is an error
				char preamble[STRSIZE];
				sprintf(preamble, "There was an error reading from %s", avg_outfile_sorted.c_str());
				perror(preamble);
				return 1;
			}
			fclose(avg_outfile_ptr);			
		
			// Print the Funseq scores to a new file that will replace the old one
			string funseq_outfile = outdir + "/permutation_" + string(perm_num) + ".funseq.txt";
			FILE *funseq_outfile_ptr = fopen(funseq_outfile.c_str(), "w");
		
			// Funseq score index
			unsigned int fs_index = 0;
		
			// Need to read the other columns from the original outfile
			// linebuf3[STRSIZE];
			outfile_ptr = fopen(outfile.c_str(), "r");
			while (fgets(linebuf3, STRSIZE, outfile_ptr) != NULL) {
				string line = string(linebuf3);
			
				// Extract chromosome, start, end (first 3 columns)
				vector<string> vec;
				for (int i = 0; i < 3; i++) {
					size_t ws_index = line.find_first_of("\t\n");
					string in = line.substr(0, ws_index);
					vec.push_back(in);
					line = line.substr(ws_index+1);
				}
			
				fprintf(funseq_outfile_ptr, "%s\t%s\t%s\t%e\n", vec[0].c_str(), vec[1].c_str(), vec[2].c_str(), funseq_scores[fs_index]);
				fs_index++;
			}
			// Check feof of vfile
			if (feof(outfile_ptr)) { // We're good
				fclose(outfile_ptr);
			} else { // It's an error
				char errstring[STRSIZE];
				sprintf(errstring, "Error reading from %s", outfile.c_str());
				perror(errstring);
				MPI_Abort(MPI_COMM_WORLD, 1);
				return 1;
			}

	// 			for (unsigned int k = 0; k < funseq_scores.size(); k++) {
	// 				fprintf(funseq_outfile_ptr, "%s\t%s\t%s\t%s\t%s\t%e\n", permuted_set[k][0].c_str(), permuted_set[k][1].c_str(), permuted_set[k][2].c_str(), permuted_set[k][3].c_str(), permuted_set[k][4].c_str(), funseq_scores[k]);
	// 			}
			fclose(funseq_outfile_ptr);

			string file_switch_1 = "rm " + outfile;
			system(file_switch_1.c_str());
			string file_switch_2 = "mv " + funseq_outfile + " " + outfile;
			system(file_switch_2.c_str());

			// DEBUG
			// printf("Filesystem magic done\n");

			// Clean up temporary folder
			string rm_com = "rm -rf " + funseq_outdir;
			system(rm_com.c_str());
		}
		
		// DEBUG
		printf("Breakpoint 5\n");
			
	} else { // Child process, do a subdivision of the work
	
		srand(0);
		
		// Receive the trimer boolean flag
		int trimer;
		MPI_Bcast(&trimer, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
		// Receive directory with wg FASTA files
		char fasta_dir_cstr[STRSIZE];
		int strlen;
		MPI_Bcast(&strlen, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(fasta_dir_cstr, strlen, MPI_CHAR, 0, MPI_COMM_WORLD);
		string fasta_dir = string(fasta_dir_cstr);
		// string(argv[5]);
		
		// Flag that indicates if all permutations are complete
		// 1 = permutations to do, 0 = all permutations done
		// Until all permutations done, loop
		int permutation_flag;
		MPI_Bcast(&permutation_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		while (permutation_flag) {
		
			// The output vector
			vector<vector<string> > permuted_set;
			
			// Vector of var_array indices for sample IDs
			vector<unsigned int> id_array;
			// vector<unsigned int> id_array_2;
		
			// Broadcast that I'm available to work
			int available_flag = 1;
			MPI_Send(&available_flag, 1, MPI_INT, 0, 9, MPI_COMM_WORLD);
		
			// Listen for the flag to receive or send
			// 0 = child receive, 1 = child send
			int flag;
			MPI_Recv(&flag, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);
	
			// BEGIN CHILD CODE
			while (flag == 0) {
		
				unsigned int variant_pointer = 0;
		
				// Begin receiving chromosome annotation data
		
				// Chromosome number
				int chr_num;
				MPI_Recv(&chr_num, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			
				// Annotation coordinates
				int chr_ann_coor_size;
				MPI_Recv(&chr_ann_coor_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
			
				int *chr_ann_coor = (int *)malloc(chr_ann_coor_size*sizeof(int));
				MPI_Recv(chr_ann_coor, chr_ann_coor_size, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
			
				// Variant coordinates
				int chr_var_coor_size;
				MPI_Recv(&chr_var_coor_size, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
			
				int *chr_var_coor = (int *)malloc(chr_var_coor_size*sizeof(int));
				MPI_Recv(chr_var_coor, chr_var_coor_size, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);
				
// 				char *chr_var_alleles = (char *)malloc(2*chr_var_coor_size*(sizeof(char)));
// 				MPI_Recv(chr_var_alleles, 2*chr_var_coor_size, MPI_CHAR, 0, 20, MPI_COMM_WORLD, &status);
			
				// Turn into the expected data structures (vectors)
				string chr = int2chr(chr_num);
				
				// DEBUG
				printf("Breakpoint Alpha\n");
			
				vector<vector<string> > var_array;
			
				for (int i = 0; i < chr_var_coor_size; i++) {
					char start_str[STRSIZE];
					sprintf(start_str, "%d", chr_var_coor[i]-1);
				
					char end_str[STRSIZE];
					sprintf(end_str, "%d", chr_var_coor[i]);
					
// 					char ref_str[STRSIZE];
// 					sprintf(ref_str, "%c", chr_var_alleles[2*i]);
// 					
// 					char alt_str[STRSIZE];
// 					sprintf(alt_str, "%c", chr_var_alleles[2*i+1]);
				
					vector<string> vec;
					vec.push_back(chr);
					vec.push_back(string(start_str));
					vec.push_back(string(end_str));
// 					vec.push_back(string(ref_str));
// 					vec.push_back(string(alt_str));
					var_array.push_back(vec);
				}
			
				vector<vector<string> > ann_array;
			
				for (int i = 0; i < chr_ann_coor_size/2; i++) {
					char start_str[STRSIZE];
					sprintf(start_str, "%d", chr_ann_coor[2*i]);
				
					char end_str[STRSIZE];
					sprintf(end_str, "%d", chr_ann_coor[2*i+1]);
				
					vector<string> vec;
					vec.push_back(chr);
					vec.push_back(string(start_str));
					vec.push_back(string(end_str));
					ann_array.push_back(vec);
				}
				
				// DEBUG
				printf("Breakpoint Beta\n");
				
				string chr_nt = "";
			
				// if (trimer) {
					// Load FASTA
					string filename = fasta_dir + "/" + chr + ".fa";
					FILE *fasta_ptr = fopen(filename.c_str(), "r");
		
					int first = 1;
					chr_nt = "";
					char linebuf_cstr[STRSIZE];
					while (fgets(linebuf_cstr, STRSIZE, fasta_ptr) != NULL) {
						string linebuf = string(linebuf_cstr);
						linebuf.erase(linebuf.find_last_not_of(" \n\r\t")+1);
						if (first) {
							first = 0;
							continue;
						}
						chr_nt += linebuf;
					}
					// Check feof of fasta_ptr
					if (feof(fasta_ptr)) { // We're good
						fclose(fasta_ptr);
					} else { // It's an error
						char errstring[STRSIZE];
						sprintf(errstring, "Error reading from %s", filename.c_str());
						perror(errstring);
						MPI_Abort(MPI_COMM_WORLD, 1);
						return 1;
					}
				// }
				
				// DEBUG
				printf("Breakpoint Ceti\n");
			
				for (unsigned int j = 0; j < ann_array.size(); j++) {
					int rand_range_start = atoi(ann_array[j][1].c_str());
					int rand_range_end = atoi(ann_array[j][2].c_str());
	
					vector<string> rand_range;
					rand_range.push_back(ann_array[j][0]);
				
					char rand_range_start_cstr[STRSIZE];
					sprintf(rand_range_start_cstr, "%d", rand_range_start);
					rand_range.push_back(string(rand_range_start_cstr));
	
					char rand_range_end_cstr[STRSIZE];
					sprintf(rand_range_end_cstr, "%d", rand_range_end);
					rand_range.push_back(string(rand_range_end_cstr));
				
					pair<unsigned int,unsigned int> range = intersecting_variantsV(var_array, rand_range, variant_pointer);
					variant_pointer = range.first;
					
					// Store the var_array indices of the variants that we're working with
					// Used to join with extra columns in output step
					vector<unsigned int> temp_id_array;
					
					// Add to id_array
					for (unsigned int k = range.first; k < range.second; k++) {
						temp_id_array.push_back(k);
					}
				
					int var_subset_count = range.second - range.first;
					if (var_subset_count == 0) {
						continue;
					}
					
					// BEGIN 3MER CODE
					map<string,vector<int> > local_nt;
				
					if (trimer) {
					
						vector<char> base; // No treble
						base.push_back('A');
						base.push_back('G');
						base.push_back('C');
						base.push_back('T');
						
						if (trimer == 3) {
							// Gather up the locations of all confidently mapped trinucleotides (capital letters)
							// Coordinates are for the second letter in the trinucleotide (where the actual mutation is located)
							// map<string,vector<int> > local_nt;
			
							for (int z = 0; z < 4; z++) {
								for (int y = 0; y < 4; y++) {
									for (int x = 0; x < 4; x++) {
										stringstream ss;
										string cur_nt;
										ss << base[z];
										ss << base[y];
										ss << base[x];
										ss >> cur_nt;
										vector<int> temp;
										local_nt[cur_nt] = temp;
									}
								}
							}
						} else if (trimer == 5) {
							// Gather up the locations of all confidently mapped pentanucleotides (capital letters)
							// Coordinates are for the third letter in the pentanucleotide (where the actual mutation is located)
							// map<string,vector<int> > local_nt;
			
							for (int z = 0; z < 4; z++) {
								for (int y = 0; y < 4; y++) {
									for (int x = 0; x < 4; x++) {
										for (int a = 0; a < 4; a++) {
											for (int b = 0; b < 4; b++) {
												stringstream ss;
												string cur_nt;
												ss << base[z];
												ss << base[y];
												ss << base[x];
												ss << base[a];
												ss << base[b];
												ss >> cur_nt;
												vector<int> temp;
												local_nt[cur_nt] = temp;
											}
										}
									}
								}
							}
						}
	
						for (int k = rand_range_start; k < rand_range_end; k++) { // 1-based index
							char nt0, nt4;
							if (trimer == 5) {
								nt0 = toupper(chr_nt[k-3]); // 0-based index
								nt4 = toupper(chr_nt[k+1]); // 0-based index
							}
							char nt1 = toupper(chr_nt[k-2]); // 0-based index
							char nt2 = toupper(chr_nt[k-1]); // 0-based index
							char nt3 = toupper(chr_nt[k]); // 0-based index
						
							// Verify there are no invalid characters
							if (nt2 != 'A' && nt2 != 'C' && nt2 != 'G' && nt2 != 'T' && nt2 != 'N') {
								char errstring[STRSIZE];
								sprintf(errstring, "Error: Invalid character detected in FASTA file: %c. Must be one of [AGCTN].\n", nt2);
								fprintf(stderr, errstring);
								return 1;
							}
						
							stringstream ss;
							string nt;
							if (trimer == 5) {
								ss << nt0;
							}
							ss << nt1;
							ss << nt2;
							ss << nt3;
							if (trimer == 5) {
								ss << nt4;
							}
							ss >> nt;
						
							// DEBUG
	// 						printf("char1: %c\n", nt1);
	// 						printf("char2: %c\n", nt2);
	// 						printf("char3: %c\n", nt3);
	// 						printf("full string: %s\n\n", nt.c_str());
						
							local_nt[nt].push_back(k-1);
						}
					}
					
					// END 3MER CODE
					
					// DEBUG
// 					printf("Size check 1: TTT: %d\n", (int)local_nt["TTT"].size());
// 					printf("Size check 1: GTC: %d\n", (int)local_nt["GTC"].size());
// 					MPI_Abort(MPI_COMM_WORLD, 1);
// 					return 1;

				// DEBUG
				printf("Breakpoint Delta\n");
				
					// Variant processing loop
					for (unsigned int k = range.first; k < range.second; k++) {
					
						if (atoi(var_array[k][2].c_str()) > rand_range_end) {
							continue;
						}
	
						// DEBUG
						// printf("%s:%s-%s\n", var_array[k][0].c_str(), var_array[k][1].c_str(), var_array[k][2].c_str());
						int new_index;
						vector<int> pos2;
						
						// BEGIN 3MER CODE
						
						if (trimer) {
							vector<string> cur_var = var_array[k];
							char cur_nt0, cur_nt4;
							if (trimer == 5) {
								cur_nt0 = toupper(chr_nt[atoi(cur_var[2].c_str())-3]);
								cur_nt4 = toupper(chr_nt[atoi(cur_var[2].c_str())+1]);
							}
							char cur_nt1 = toupper(chr_nt[atoi(cur_var[2].c_str())-2]);
							char cur_nt2 = toupper(chr_nt[atoi(cur_var[2].c_str())-1]); // 0-based index
							char cur_nt3 = toupper(chr_nt[atoi(cur_var[2].c_str())]);
		
							stringstream ss;
							string cur_nt;
							if (trimer == 5) {
								ss << cur_nt0;
							}
							ss << cur_nt1;
							ss << cur_nt2;
							ss << cur_nt3;
							if (trimer == 5) {
								ss << cur_nt4;
							}
							ss >> cur_nt;
						
							// If there is an N in this string, we skip this variant
							if (cur_nt.find_first_of('N') != string::npos) {
								continue;
							}
		
							// DEBUG
							// printf("DEBUG: %c,%c\n", cur_nt1, cur_nt2);
							// printf("DEBUG: cur_nt: %s\n", cur_nt.c_str());
		
							vector<int> pos = local_nt[cur_nt];
						
							// DEBUG
							// printf("Size check 2: %s: %d\n", cur_nt.c_str(), (int)local_nt[cur_nt].size());
		
							// If no positions are available, end program with an error and suggest
							// a larger bin size
							if (pos.size()-1 == 0) {
								char errstring[STRSIZE];
								sprintf(errstring, "Error: No valid permutations positions for a variant in bin %s:%s-%s. Consider using a larger bin size.\n",
												ann_array[j][0].c_str(), ann_array[j][1].c_str(), ann_array[j][2].c_str());
								fprintf(stderr, errstring);
								MPI_Abort(MPI_COMM_WORLD, 1);
								return 1;
							}
		
							// vector<int> pos2;
							for (unsigned int l = 0; l < pos.size(); l++) {
								if (pos[l] != atoi(cur_var[2].c_str())-1) {
									pos2.push_back(pos[l]);
								}
							}
						
							// DEBUG
	// 						printf("Size check 3: pos: %d\n", (int)pos.size());
	// 						printf("Size check 3: pos2: %d\n", (int)pos2.size());
						
							// Pick new position
							new_index = rand() % (pos2.size()); // Selection in interval [0,pos2.size()-1]
						} else {
							do {
								new_index = rand() % (rand_range_end-rand_range_start); // Selection in interval [0,(rand_range_end-rand_range_start)-1]
							} while (chr_nt[new_index] == 'N');
						}
						
						// END 3MER CODE
						
						// DEBUG
						printf("Breakpoint Epsilon\n");
						
						vector<string> vec;
						vec.push_back(var_array[k][0]);
		
						if (trimer) {
						
							char start_cstr[STRSIZE];
							sprintf(start_cstr, "%d", pos2[new_index]); // 0-based
							vec.push_back(string(start_cstr));
		
							char end_cstr[STRSIZE];
							sprintf(end_cstr, "%d", pos2[new_index]+1); // 1-based
							vec.push_back(string(end_cstr));
							
// 							vec.push_back(var_array[k][3]);
// 							vec.push_back(var_array[k][4]);
							
						} else {
							
							char start_cstr[STRSIZE];
							sprintf(start_cstr, "%d", new_index); // 0-based
							vec.push_back(string(start_cstr));
				
							char end_cstr[STRSIZE];
							sprintf(end_cstr, "%d", new_index+1); // 1-based
							vec.push_back(string(end_cstr));
							
// 							bool is_purine; // Otherwise, pyrimidine
// 							bool is_transition; // Otherwise, transversion
// 							
// 							char old_ref = toupper(var_array[k][3][0]);
// 							char old_alt = toupper(var_array[k][4][0]);
// 							
// 							if (old_ref == 'A' || old_ref == 'G') {
// 								is_purine = true;
// 							} else {
// 								is_purine = false;
// 							}
// 					
// 							if (is_purine) {
// 								if (old_alt == 'A' || old_alt == 'G') {
// 									is_transition = true;
// 								} else {
// 									is_transition = false;
// 								}
// 							} else {
// 								if (old_alt == 'A' || old_alt == 'G') {
// 									is_transition = false;
// 								} else {
// 									is_transition = true;
// 								}
// 							}
// 							
// 							char ref = toupper(chr_nt[new_index]);
// 							string alt;
// 					
// 							if (is_transition) {
// 								if (ref == 'A') {
// 									alt = "G";
// 								} else if (ref == 'G') {
// 									alt = "A";
// 								} else if (ref == 'C') {
// 									alt = "T";
// 								} else if (ref == 'T') {
// 									alt = "C";
// 								}
// 							} else {
// 								int rando = rand() % 2;
// 								if (ref == 'A' || ref == 'G') {
// 									// Choose between C and T
// 									if (rando) {
// 										alt = "C";
// 									} else {
// 										alt = "T";
// 									}
// 								} else {
// 									// Choose between A and G
// 									if (rando) {
// 										alt = "A";
// 									} else {
// 										alt = "G";
// 									}
// 								}
// 							}
// 							
// 							char ref_str[STRSIZE];
// 							sprintf(ref_str, "%c", ref);
// 					
// 							vec.push_back(string(ref_str));
// 							vec.push_back(alt);
							
						}
		
						permuted_set.push_back(vec);
						id_array.push_back(temp_id_array[k]);
					}
				}
				free(chr_ann_coor);
				free(chr_var_coor);
				// free(chr_var_alleles);
				
				MPI_Send(&available_flag, 1, MPI_INT, 0, 9, MPI_COMM_WORLD);
				MPI_Recv(&flag, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);
			}
			
			// DEBUG
			printf("Breakpoint Gamma\n");
		
			// Sending time
			// Send chr and end coordinates of the permuted variants
			int permuted_var_coor_size = 2*permuted_set.size();
			MPI_Send(&permuted_var_coor_size, 1, MPI_INT, 0, 6, MPI_COMM_WORLD);
			
			// DEBUG
			printf("permuted_var_coor_size: %d\n", permuted_var_coor_size);
			
			// Proceed if nonzero size
			if (permuted_var_coor_size > 0) {
				int *permuted_var_coor = (int *)malloc(permuted_var_coor_size*sizeof(int));
				// char *permuted_var_alleles = (char *)malloc(permuted_var_coor_size*sizeof(char));
				unsigned int *permuted_var_id = (unsigned int *)malloc((permuted_var_coor_size/2)*sizeof(unsigned int));
		
				for (unsigned int i = 0; i < permuted_set.size(); i++) {
					permuted_var_coor[2*i] = chr2int(permuted_set[i][0]);
					permuted_var_coor[2*i+1] = atoi(permuted_set[i][2].c_str());
					// permuted_var_alleles[2*i] = permuted_set[i][3][0];
					// permuted_var_alleles[2*i+1] = permuted_set[i][4][0];
					permuted_var_id[i] = id_array[i];
				}
				
				// DEBUG
				printf("permuted_var_coor[0]: %d\n", permuted_var_coor[0]);
				printf("permuted_var_coor[1]: %d\n", permuted_var_coor[1]);
				printf("permuted_var_id[0]: %d\n", permuted_var_id[0]);
				printf("permuted_var_id[1]: %d\n", permuted_var_id[1]);
				printf("permuted_var_id[2]: %d\n", permuted_var_id[2]);
		
				MPI_Send(permuted_var_coor, permuted_var_coor_size, MPI_INT, 0, 7, MPI_COMM_WORLD);
				// MPI_Send(permuted_var_alleles, permuted_var_coor_size, MPI_CHAR, 0, 21, MPI_COMM_WORLD);
				MPI_Send(permuted_var_id, (permuted_var_coor_size/2), MPI_UNSIGNED, 0, 23, MPI_COMM_WORLD);
				
				// DEBUG
				printf("Breakpoint Gamma Quadrant\n");
			
				free(permuted_var_coor);
				// free(permuted_var_alleles);
				free(permuted_var_id);
			}
			
			MPI_Bcast(&permutation_flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
			// END CHILD CODE
			
			// DEBUG
			printf("Breakpoint Hydra\n");
		}
	}

	// DEBUG
	printf("Breakpoint 3\n");
	
	// Verdun
	MPI_Finalize();
	return 0;
}
