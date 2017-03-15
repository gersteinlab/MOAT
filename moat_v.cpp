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
#include "variant_permutation_v3.h"

using namespace std;

#define STRSIZE 256
#define BIGSTRSIZE 1024

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
	
	/* User-supplied arguments */
	
	// Is the trinucleotide preservation option enabled?
	bool trimer;
	
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
	
	// Option to specify whether to calculate Funseq scores on the permutations
	bool funseq_opt;
	
	if (argc != 10) {
		fprintf(stderr, "Usage: moat_v [3mer preservation option (y/n)] [# permuted datasets] [permutation window radius] [min width] [prohibited regions file] [FASTA dir] [variant file] [output directory] [Funseq option (y/n)]. Exiting.\n");
		return 1;
	} else {
	
		if (argv[1][0] == 'y') {
			trimer = true;
		} else if (argv[1][0] == 'n') {
			trimer = false;
		} else {
			fprintf(stderr, "Invalid option for 3mer preservation option: \'%c\'. Must be either \'y\' or \'n\'. Exiting.\n", argv[1][0]);
			return 1;
		}
			
		num_permutations = atoi(argv[2]);
		window_radius = atoi(argv[3]);
		min_width = atoi(argv[4]);
		prohibited_file = string(argv[5]);
		fasta_dir = string(argv[6]);
		vfile = string(argv[7]);
		outdir = string(argv[8]);
		
		if (argv[9][0] == 'y') {
			funseq_opt = true;
		} else if (argv[9][0] == 'n') {
			funseq_opt = false;
		} else {
			fprintf(stderr, "Invalid option for Funseq option: \'%c\'. Must be either \'y\' or \'n\'. Exiting.\n", argv[9][0]);
			return 1;
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
	if (pbuf.st_size == 0) {
		fprintf(stderr, "Error: Prohibited regions file cannot be empty. Exiting.\n");
		return 1;
	}
	
	// Check that the FASTA directory is a valid path
	struct stat fbuf;
	if (stat(fasta_dir.c_str(), &fbuf)) { // Report the error and exit
		fprintf(stderr, "Error trying to stat %s: %s\n", fasta_dir.c_str(), strerror(errno));
		return 1;
	}
	
	// Check that the outdir is a valid path
	struct stat obuf;
	if (stat(outdir.c_str(), &obuf)) { // Report the error and exit
		fprintf(stderr, "Error trying to stat %s: %s\n", outdir.c_str(), strerror(errno));
		return 1;
	}
	
	// Check that Funseq2 is installed, if it is to be used
	if (funseq_opt) {
		string test_funseq2 = "command -v funseq2.sh > /dev/null";
		if (system(test_funseq2.c_str())) {
			fprintf(stderr, "Error: Funseq2 is not installed. Please install Funseq2 ");
			fprintf(stderr, "(funseq2.gersteinlab.org) and ensure the install path has been ");
			fprintf(stderr, "added to your PATH environment variable.\n");
			return 1;
		}
	}
	
	/* Data structures for the starting data */
	// Variant array, contains variants of the format vector(chr, start, end, ref allele, alt allele)
	vector<vector<string> > var_array;
	
	// Annotation array, contains annotations of the format vector(chr, start, end, ann_name)
	// Will be generated from the whole genome coordinates at runtime
	vector<vector<string> > ann_array;
	
	// Prohibited regions array, contains annotations of the format vector(chr, start, end)
	vector<vector<string> > prohibited_regions;
	
	// Bring variant file data into memory
	// Save the first 5 columns, ignore the rest if there are any
	char linebuf[STRSIZE];
	FILE *vfile_ptr = fopen(vfile.c_str(), "r");
	while (fgets(linebuf, STRSIZE, vfile_ptr) != NULL) {
		string line = string(linebuf);
		
		// DEBUG
		// printf("%s\n", line.c_str());
		
		// Extract chromosome, start, end, ref and alt from line (first 5 columns)
		vector<string> vec;
		for (int i = 0; i < 5; i++) {
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
	// printf("Breakpoint 1\n");
	
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
// 	FILE *testfile_ptr = fopen("test-bin-code/testfile.txt", "w");
// 	for (unsigned int i = 0; i < ann_array.size(); i++) {
// 		fprintf(testfile_ptr, "%s\t%s\t%s\n", ann_array[i][0].c_str(), ann_array[i][1].c_str(), ann_array[i][2].c_str());
// 	}
// 	fclose(testfile_ptr);
// 	return 0;
	
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
	// printf("Breakpoint 2\n");
	
	// DEBUG - check ann_array values after prohibited region subtraction
// 	FILE *testfile_ptr = fopen("test-bin-code/testfile.txt", "w");
// 	for (unsigned int i = 0; i < ann_array.size(); i++) {
// 		fprintf(testfile_ptr, "%s\t%s\t%s\n", ann_array[i][0].c_str(), ann_array[i][1].c_str(), ann_array[i][2].c_str());
// 	}
// 	fclose(testfile_ptr);
// 	return 0;
	
	FILE *fasta_ptr = NULL;
	string last_chr = "";
	// int char_pointer;
	string chr_nt;
	
	/* Permutate variant locations */
	srand(0);
	
	for (int i = 0; i < num_permutations; i++) {
		
		// Open new output file
		char perm_num[STRSIZE];
		sprintf(perm_num, "%d", i+1);
		
		string outfile = outdir + "/permutation_" + string(perm_num) + ".txt";
		FILE *outfile_ptr = fopen(outfile.c_str(), "w");
	
		unsigned int variant_pointer = 0;
		// vector<int> this_permutation_counts;
		
		for (unsigned int j = 0; j < ann_array.size(); j++) {
		
			vector<vector<string> > permuted_set;
		
			if (trimer) {
				if (last_chr != ann_array[j][0]) {
					// newFastaFilehandle(fasta_ptr, ann_array[j][0]);
					// char_pointer = 0;
				
					string filename = fasta_dir + "/" + ann_array[j][0] + ".fa";
					fasta_ptr = fopen(filename.c_str(), "r");
				
					int first = 1;
					last_chr = ann_array[j][0];
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
						return 1;
					}
				}
			}
			
			// DEBUG
// 			printf("Char 11998: %c\n", chr_nt[11997]);
// 			printf("Char 11999: %c\n", chr_nt[11998]);
// 			printf("Char 12000: %c\n", chr_nt[11999]);
// 			printf("Char 12001: %c\n", chr_nt[12000]);
// 			printf("Char 12002: %c\n", chr_nt[12001]);
			
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
			
			pair<unsigned int,unsigned int> range = intersecting_variants(var_array, rand_range, variant_pointer);
			variant_pointer = range.first;
			
			// DEBUG
			// printf("Breakpoint 3\n");
			// printf("%d,%d\n", range.first, range.second);
			
			int var_subset_count = range.second - range.first;
			if (var_subset_count == 0) {
				continue;
			}
			
			// BEGIN 3MER CODE
			map<string,vector<int> > local_nt;
			
			if (trimer) {
			
				// Gather up the locations of all confidently mapped trinucleotides (capital letters)
				// Coordinates are for the second letter in the trinucleotide (where the actual mutation is located)
				vector<int> AAA;
				vector<int> AAG;
				vector<int> AAC;
				vector<int> AAT;
			
				vector<int> AGA;
				vector<int> AGG;
				vector<int> AGC;
				vector<int> AGT;
			
				vector<int> ACA;
				vector<int> ACG;
				vector<int> ACC;
				vector<int> ACT; // Download from act.gersteinlab.org today!
			
				vector<int> ATA;
				vector<int> ATG;
				vector<int> ATC;
				vector<int> ATT; // Verizon's not going to be happy about this
			
				vector<int> GAA;
				vector<int> GAG;
				vector<int> GAC;
				vector<int> GAT;
			
				vector<int> GGA;
				vector<int> GGG; // Good game Greg
				vector<int> GGC;
				vector<int> GGT;
			
				vector<int> GCA;
				vector<int> GCG;
				vector<int> GCC; // Your friendly neighborhood C compiler
				vector<int> GCT;
			
				vector<int> GTA;
				vector<int> GTG;
				vector<int> GTC;
				vector<int> GTT;
			
				vector<int> CAA;
				vector<int> CAG; // Report to Lee Adama for mission briefing
				vector<int> CAC;
				vector<int> CAT; // Cats and dogs living together... total anarchy
			
				vector<int> CGA;
				vector<int> CGG;
				vector<int> CGC;
				vector<int> CGT;
			
				vector<int> CCA;
				vector<int> CCG;
				vector<int> CCC;
				vector<int> CCT;
			
				vector<int> CTA;
				vector<int> CTG;
				vector<int> CTC;
				vector<int> CTT;
			
				vector<int> TAA;
				vector<int> TAG; // You're it!
				vector<int> TAC;
				vector<int> TAT;
			
				vector<int> TGA;
				vector<int> TGG;
				vector<int> TGC;
				vector<int> TGT;
			
				vector<int> TCA; // Is your favorite show on the bubble?
				vector<int> TCG; // With one more letter, it could've been MTG
				vector<int> TCC;
				vector<int> TCT;
			
				vector<int> TTA;
				vector<int> TTG;
				vector<int> TTC; // The words of the prophet are written on the TTC
				vector<int> TTT;
			
				local_nt["AAA"] = AAA;
				local_nt["AAG"] = AAG;
				local_nt["AAC"] = AAC;
				local_nt["AAT"] = AAT;
			
				local_nt["AGA"] = AGA;
				local_nt["AGG"] = AGG;
				local_nt["AGC"] = AGC;
				local_nt["AGT"] = AGT;
			
				local_nt["ACA"] = ACA;
				local_nt["ACG"] = ACG;
				local_nt["ACC"] = ACC;
				local_nt["ACT"] = ACT;
			
				local_nt["ATA"] = ATA;
				local_nt["ATG"] = ATG;
				local_nt["ATC"] = ATC;
				local_nt["ATT"] = ATT;
			
				local_nt["GAA"] = GAA;
				local_nt["GAG"] = GAG;
				local_nt["GAC"] = GAC;
				local_nt["GAT"] = GAT;
			
				local_nt["GGA"] = GGA;
				local_nt["GGG"] = GGG;
				local_nt["GGC"] = GGC;
				local_nt["GGT"] = GGT;
			
				local_nt["GCA"] = GCA;
				local_nt["GCG"] = GCG;
				local_nt["GCC"] = GCC;
				local_nt["GCT"] = GCT;
			
				local_nt["GTA"] = GTA;
				local_nt["GTG"] = GTG;
				local_nt["GTC"] = GTC;
				local_nt["GTT"] = GTT;
			
				local_nt["CAA"] = CAA;
				local_nt["CAG"] = CAG;
				local_nt["CAC"] = CAC;
				local_nt["CAT"] = CAT;
			
				local_nt["CGA"] = CGA;
				local_nt["CGG"] = CGG;
				local_nt["CGC"] = CGC;
				local_nt["CGT"] = CGT;
			
				local_nt["CCA"] = CCA;
				local_nt["CCG"] = CCG;
				local_nt["CCC"] = CCC;
				local_nt["CCT"] = CCT;
			
				local_nt["CTA"] = CTA;
				local_nt["CTG"] = CTG;
				local_nt["CTC"] = CTC;
				local_nt["CTT"] = CTT;
			
				local_nt["TAA"] = TAA;
				local_nt["TAG"] = TAG;
				local_nt["TAC"] = TAC;
				local_nt["TAT"] = TAT;
			
				local_nt["TGA"] = TGA;
				local_nt["TGG"] = TGG;
				local_nt["TGC"] = TGC;
				local_nt["TGT"] = TGT;
			
				local_nt["TCA"] = TCA;
				local_nt["TCG"] = TCG;
				local_nt["TCC"] = TCC;
				local_nt["TCT"] = TCT;
			
				local_nt["TTA"] = TTA;
				local_nt["TTG"] = TTG;
				local_nt["TTC"] = TTC;
				local_nt["TTT"] = TTT;
			
				for (int k = rand_range_start; k < rand_range_end; k++) { // 1-based index
			
					// Don't read in characters if it will read off either end
					if (k == 1 || k == hg19_coor[ann_array[j][0]]) {
						continue;
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
				
					if (nt1 == 'A' && nt2 == 'A' && nt3 == 'A') {
						local_nt["AAA"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'A' && nt3 == 'G') {
						local_nt["AAG"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'A' && nt3 == 'C') {
						local_nt["AAC"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'A' && nt3 == 'T') {
						local_nt["AAT"].push_back(k-1);
					
					} else if (nt1 == 'A' && nt2 == 'G' && nt3 == 'A') {
						local_nt["AGA"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'G' && nt3 == 'G') {
						local_nt["AGG"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'G' && nt3 == 'C') {
						local_nt["AGC"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'G' && nt3 == 'T') {
						local_nt["AGT"].push_back(k-1);
					
					} else if (nt1 == 'A' && nt2 == 'C' && nt3 == 'A') {
						local_nt["ACA"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'C' && nt3 == 'G') {
						local_nt["ACG"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'C' && nt3 == 'C') {
						local_nt["ACC"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'C' && nt3 == 'T') {
						local_nt["ACT"].push_back(k-1);
					
					} else if (nt1 == 'A' && nt2 == 'T' && nt3 == 'A') {
						local_nt["ATA"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'T' && nt3 == 'G') {
						local_nt["ATG"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'T' && nt3 == 'C') {
						local_nt["ATC"].push_back(k-1);
					} else if (nt1 == 'A' && nt2 == 'T' && nt3 == 'T') {
						local_nt["ATT"].push_back(k-1);
					
					} else if (nt1 == 'G' && nt2 == 'A' && nt3 == 'A') {
						local_nt["GAA"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'A' && nt3 == 'G') {
						local_nt["GAG"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'A' && nt3 == 'C') {
						local_nt["GAC"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'A' && nt3 == 'T') {
						local_nt["GAT"].push_back(k-1);
					
					} else if (nt1 == 'G' && nt2 == 'G' && nt3 == 'A') {
						local_nt["GGA"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'G' && nt3 == 'G') {
						local_nt["GGG"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'G' && nt3 == 'C') {
						local_nt["GGC"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'G' && nt3 == 'T') {
						local_nt["GGT"].push_back(k-1);
					
					} else if (nt1 == 'G' && nt2 == 'C' && nt3 == 'A') {
						local_nt["GCA"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'C' && nt3 == 'G') {
						local_nt["GCG"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'C' && nt3 == 'C') {
						local_nt["GCC"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'C' && nt3 == 'T') {
						local_nt["GCT"].push_back(k-1);
					
					} else if (nt1 == 'G' && nt2 == 'T' && nt3 == 'A') {
						local_nt["GTA"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'T' && nt3 == 'G') {
						local_nt["GTG"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'T' && nt3 == 'C') {
						local_nt["GTC"].push_back(k-1);
					} else if (nt1 == 'G' && nt2 == 'T' && nt3 == 'T') {
						local_nt["GTT"].push_back(k-1);
				
					} else if (nt1 == 'C' && nt2 == 'A' && nt3 == 'A') {
						local_nt["CAA"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'A' && nt3 == 'G') {
						local_nt["CAG"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'A' && nt3 == 'C') {
						local_nt["CAC"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'A' && nt3 == 'T') {
						local_nt["CAT"].push_back(k-1);
					
					} else if (nt1 == 'C' && nt2 == 'G' && nt3 == 'A') {
						local_nt["CGA"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'G' && nt3 == 'G') {
						local_nt["CGG"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'G' && nt3 == 'C') {
						local_nt["CGC"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'G' && nt3 == 'T') {
						local_nt["CGT"].push_back(k-1);
					
					} else if (nt1 == 'C' && nt2 == 'C' && nt3 == 'A') {
						local_nt["CCA"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'C' && nt3 == 'G') {
						local_nt["CCG"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'C' && nt3 == 'C') {
						local_nt["CCC"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'C' && nt3 == 'T') {
						local_nt["CCT"].push_back(k-1);
					
					} else if (nt1 == 'C' && nt2 == 'T' && nt3 == 'A') {
						local_nt["CTA"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'T' && nt3 == 'G') {
						local_nt["CTG"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'T' && nt3 == 'C') {
						local_nt["CTC"].push_back(k-1);
					} else if (nt1 == 'C' && nt2 == 'T' && nt3 == 'T') {
						local_nt["CTT"].push_back(k-1);
					
					} else if (nt1 == 'T' && nt2 == 'A' && nt3 == 'A') {
						local_nt["TAA"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'A' && nt3 == 'G') {
						local_nt["TAG"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'A' && nt3 == 'C') {
						local_nt["TAC"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'A' && nt3 == 'T') {
						local_nt["TAT"].push_back(k-1);
					
					} else if (nt1 == 'T' && nt2 == 'G' && nt3 == 'A') {
						local_nt["TGA"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'G' && nt3 == 'G') {
						local_nt["TGG"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'G' && nt3 == 'C') {
						local_nt["TGC"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'G' && nt3 == 'T') {
						local_nt["TGT"].push_back(k-1);
					
					} else if (nt1 == 'T' && nt2 == 'C' && nt3 == 'A') {
						local_nt["TCA"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'C' && nt3 == 'G') {
						local_nt["TCG"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'C' && nt3 == 'C') {
						local_nt["TCC"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'C' && nt3 == 'T') {
						local_nt["TCT"].push_back(k-1);
					
					} else if (nt1 == 'T' && nt2 == 'T' && nt3 == 'A') {
						local_nt["TTA"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'T' && nt3 == 'G') {
						local_nt["TTG"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'T' && nt3 == 'C') {
						local_nt["TTC"].push_back(k-1);
					} else if (nt1 == 'T' && nt2 == 'T' && nt3 == 'T') {
						local_nt["TTT"].push_back(k-1);
					}
				}
			}
			
			// END 3MER CODE
			
			// Variant processing loop
			for (unsigned int k = range.first; k < range.second; k++) {
			
				// DEBUG
				// printf("%s:%s-%s\n", var_array[k][0].c_str(), var_array[k][1].c_str(), var_array[k][2].c_str());
				
				int new_index;
				vector<int> pos2;
				
				// BEGIN 3MER CODE
				if (trimer) {
			
					vector<string> cur_var = var_array[k];
					char cur_nt1 = toupper(chr_nt[atoi(cur_var[2].c_str())-2]);
					char cur_nt2 = toupper(chr_nt[atoi(cur_var[2].c_str())-1]); // 0-based index
					char cur_nt3 = toupper(chr_nt[atoi(cur_var[2].c_str())]);
				
					stringstream ss;
					string cur_nt;
					ss << cur_nt1;
					ss << cur_nt2;
					ss << cur_nt3;
					ss >> cur_nt;
				
					// If there is an N in this string, we skip this variant
					if (cur_nt.find_first_of('N') != string::npos) {
						continue;
					}
				
					// DEBUG
					// printf("DEBUG: %c,%c\n", cur_nt1, cur_nt2);
					// printf("DEBUG: cur_nt: %s\n", cur_nt.c_str());
				
					vector<int> pos = local_nt[cur_nt];
				
					// If no positions are available, end program with an error and suggest
					// a larger bin size
					if (pos.size()-1 == 0) {
						char errstring[STRSIZE];
						sprintf(errstring, "Error: No valid permutations positions for a variant in bin %s:%s-%s. Consider using a larger bin size.\n",
										ann_array[j][0].c_str(), ann_array[j][1].c_str(), ann_array[j][2].c_str());
						fprintf(stderr, errstring);
						return 1;
					}
				
					// vector<int> pos2;
					for (unsigned int l = 0; l < pos.size(); l++) {
						if (pos[l] != atoi(cur_var[2].c_str())-1) {
							pos2.push_back(pos[l]);
						}
					}
					// Pick new position
					new_index = rand() % (pos2.size()); // Selection in interval [0,pos2.size()-1]
				} else {
					new_index = rand() % (rand_range_end-rand_range_start); // Selection in interval [0,(rand_range_end-rand_range_start)-1]
				}
				
				// END 3MER CODE
				
				vector<string> vec;
				vec.push_back(var_array[k][0]);
				
				if (trimer) {
				
					char start_cstr[STRSIZE];
					sprintf(start_cstr, "%d", pos2[new_index]); // 0-based
					vec.push_back(string(start_cstr));
				
					char end_cstr[STRSIZE];
					sprintf(end_cstr, "%d", pos2[new_index]+1); // 1-based
					vec.push_back(string(end_cstr));
					
				} else {
				
					char start_cstr[STRSIZE];
					sprintf(start_cstr, "%d", new_index); // 0-based
					vec.push_back(string(start_cstr));
				
					char end_cstr[STRSIZE];
					sprintf(end_cstr, "%d", new_index+1); // 1-based
					vec.push_back(string(end_cstr));
					
				}
				
				vec.push_back(var_array[k][3]);
				vec.push_back(var_array[k][4]);
				
				permuted_set.push_back(vec);
			}
			
			// Read in the basepairs we need for this region
// 			string local_nt = "";
// 			if (remainder != "") {
// 				local_nt += remainder;
// 			}
// 			local_nt += fillBuffer(fasta_ptr, &char_pointer, rand_range_start, rand_range_end, &remainder);
			// Need to save the trailing bit of the strbuf
			
			// END NEW CODE
			
			// vector<vector<string> > permuted_set = permute_variants(var_subset_count, rand_range);
			sort(permuted_set.begin(), permuted_set.end(), cmpIntervals);
		
			// DEBUG
			// printf("Loop iter %d\n", i);
	// 		for (unsigned int j = 0; j < permuted_set.size(); j++) {
	// 			printf("Variant from set %d: %s %d %d\n", i, (permuted_set[j][0]).c_str(), atoi((permuted_set[j][1]).c_str()), atoi((permuted_set[j][2]).c_str()));
	// 		}
	// 		return 0;
			// printf("Permutation %d, Permuted set size: %d\n", i+1, (int)permuted_set.size());
		
			for (unsigned int k = 0; k < permuted_set.size(); k++) {
				fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%s\n", permuted_set[k][0].c_str(), permuted_set[k][1].c_str(), permuted_set[k][2].c_str(), permuted_set[k][3].c_str(), permuted_set[k][4].c_str());
			}
			last_chr = ann_array[j][0];
		}
		fclose(outfile_ptr);
		
		// Data available for Funseq step
		vector<double> funseq_scores;
		if (funseq_opt) {
			
			string funseq_loc = exec("command -v funseq2.sh");
			size_t index = funseq_loc.find_last_of("/");
			funseq_loc = funseq_loc.substr(0, index);
		
			// Retrieve current working directory for temporary Funseq2 output
			string funseq_outdir = exec("pwd");
			funseq_outdir.erase(funseq_outdir.find_last_not_of(" \n\r\t")+1);
			funseq_outdir += "/funseq";
			
			// Verify that funseq output directory exists, or create it if it doesn't
			string command = "mkdir -p " + funseq_outdir;
			system(command.c_str());
			
			// Convert outfile to absolute path, if it is not already
			char abs_path_outfile[PATH_MAX];
			errno = 0;
			realpath(outfile.c_str(), abs_path_outfile);
			if (errno) {
				fprintf(stderr, "Error resolving absolute path of permutation variant file: %s\n", strerror(errno));
				fprintf(stderr, "Exiting.\n");
				return 1;
			}
		
			string funseq2_command = "cd " + funseq_loc + "; ./funseq2.sh -f " + abs_path_outfile + " -inf bed -outf bed -o " + funseq_outdir;
			system(funseq2_command.c_str());
			
			// Collect sum of Funseq scores per annotation
			vector<vector<string> > funseq_output;
		
			// Read in "Output.bed"
			int first = 1;
			char linebuf2[BIGSTRSIZE];
			string funseq_output_file = funseq_outdir + "/Output.bed";
			FILE *ffile_ptr = fopen(funseq_output_file.c_str(), "r");
			while (fgets(linebuf2, BIGSTRSIZE, ffile_ptr) != NULL) {
		
				if (first) {
					first = 0;
					continue;
				}
		
				string line = string(linebuf2);
			
				vector<string> vec;
				for (int i = 0; i < 7; i++) {
					size_t ws_index = line.find_first_of("\t\n");
					string in = line.substr(0, ws_index);
					vec.push_back(in);
					line = line.substr(ws_index+1);
				}
			
				// If this is not a standard chromosome, then remove this row
				if (chr2int(vec[0]) == 0) {
					continue;
				}
			
				funseq_output.push_back(vec);
			}
			// Check feof of vfile
			if (feof(ffile_ptr)) { // We're good
				fclose(ffile_ptr);
			} else { // It's an error
				char errstring[STRSIZE];
				sprintf(errstring, "Error reading from %s", funseq_output_file.c_str());
				perror(errstring);
				return 1;
			}
		
			// Sort
			sort(funseq_output.begin(), funseq_output.end(), cmpIntervals);
		
			// Gather up and sum the Funseq values over each annotation
			// unsigned int funseq_var_pointer = 0;
			for (unsigned int j = 0; j < funseq_output.size(); j++) {
				string info_str = funseq_output[j][6];
				double funseq_sum = 0.0;
				double coding_score;
				double nc_score;
				for (int i = 0; i < 14; i++) {
					size_t ws_index = info_str.find_first_of(";");
					string in = info_str.substr(0, ws_index);
					info_str = info_str.substr(ws_index+1);
					if (i == 12) {
						if (in != ".") {
							coding_score = atof(in.c_str());
						} else {
							coding_score = -1.0;
						}
					} else if (i == 13) {
						if (in != ".") {
							nc_score = atof(in.c_str());
						} else {
							nc_score = -1.0;
						}
					}
				}
			
				if (coding_score != -1.0) {
					funseq_sum = coding_score;
				} else {
					funseq_sum = nc_score;
				}
				funseq_scores.push_back(funseq_sum);
			}
		}
			
		// Print the Funseq scores to a new file that will replace the old one
		string funseq_outfile = outdir + "/permutation_" + string(perm_num) + ".funseq.txt";
		FILE *funseq_outfile_ptr = fopen(funseq_outfile.c_str(), "w");
		
		for (unsigned int k = 0; k < funseq_scores.size(); k++) {
			fprintf(funseq_outfile_ptr, "%s\t%s\t%s\t%s\t%s\t%e\n", permuted_set[k][0].c_str(), permuted_set[k][1].c_str(), permuted_set[k][2].c_str(), permuted_set[k][3].c_str(), permuted_set[k][4].c_str(), funseq_scores[k]);
		}
		fclose(funseq_outfile_ptr);
		
		string file_switch_1 = "rm " + outfile;
		system(file_switch_1.c_str());
		string file_switch_2 = "mv " + funseq_outfile + " " + outfile;
		system(file_switch_2.c_str());
	}
	// }
	
	// DEBUG
	// printf("Breakpoint 3\n");
	
	// Verdun
	return 0;
}
