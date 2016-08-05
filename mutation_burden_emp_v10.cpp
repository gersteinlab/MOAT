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
#include <math.h>
#include <stdexcept>
#include <float.h>
#include <limits.h>
#include "variant_permutation_v3.h"

using namespace std;

#define STRSIZE 256

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

// In v10, covariate inputs have been added, and the bins are clustered on these
// criteria, and all the variants in these clusters are permuted within the cluster bins

/* 
 * Depends on bigWigAverageOverBed
 */

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

// Helper function that calculates the Euclidean distance between two vectors
// Assumes the two vectors are the same length
double euclidean(vector<double> &a, vector<double> &b) {
	double sum = 0.0;
	for (unsigned int i = 0; i < a.size(); i++) {
		sum += pow((b[i]-a[i]), 2.0);
	}
	return sqrt(sum);
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
	
	// Covariate signal files in bigWig format
	vector<string> covar_files;
	
	if (argc < 9) {
		printf("Usage: mutation_burden_emp_v10 [# permuted datasets] [permutation window radius] [min width] [prohibited regions file] [FASTA dir] [variant file] [output folder] [covariate files ...]. Exiting.\n");
		return 1;
	} else {
		num_permutations = atoi(argv[1]);
		window_radius = atoi(argv[2]);
		min_width = atoi(argv[3]);
		prohibited_file = string(argv[4]);
		fasta_dir = string(argv[5]);
		vfile = string(argv[6]);
		outdir = string(argv[7]);
		
		for (int i = 8; i < argc; i++) {
			covar_files.push_back(string(argv[i]));
		}
	}
	
	// Verify files, and import data to memory
	struct stat vbuf;
	if (stat(vfile.c_str(), &vbuf)) { // Report the error and exit
		printf("Error trying to stat %s: %s\n", vfile.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (vbuf.st_size == 0) {
		printf("Error: Variant file cannot be empty. Exiting.\n");
		return 1;
	}
	
	struct stat pbuf;
	if (stat(prohibited_file.c_str(), &pbuf)) { // Report the error and exit
		printf("Error trying to stat %s: %s\n", prohibited_file.c_str(), strerror(errno));
		return 1;
	}
	// Check that the file is not empty
	if (pbuf.st_size == 0) {
		printf("Error: Prohibited regions file cannot be empty. Exiting.\n");
		return 1;
	}
	
	// Check that the FASTA directory is a valid path
	struct stat fbuf;
	if (stat(fasta_dir.c_str(), &fbuf)) { // Report the error and exit
		printf("Error trying to stat %s: %s\n", fasta_dir.c_str(), strerror(errno));
		return 1;
	}
	
	// Check that the outdir is a valid path
	struct stat obuf;
	if (stat(outdir.c_str(), &obuf)) { // Report the error and exit
		printf("Error trying to stat %s: %s\n", outdir.c_str(), strerror(errno));
		return 1;
	}
	
	// Verify that bigWigAverageOverBed is in the same directory as this program
	struct stat avgbuf;
	char avgoverbed_cstr[] = "./bigWigAverageOverBed";
	if (stat(avgoverbed_cstr, &avgbuf)) {
		printf("Error: bigWigAverageOverBed is not in the same directory. Exiting.\n");
		return 1;
	}
	
	/* Data structures for the starting data */
	// Variant array, contains variants of the format vector(chr, start, end)
	vector<vector<string> > var_array;
	
	// Annotation array, contains annotations of the format vector(chr, start, end)
	// Will be generated from the whole genome coordinates at runtime
	vector<vector<string> > ann_array;
	
	// Prohibited regions array, contains annotations of the format vector(chr, start, end)
	vector<vector<string> > prohibited_regions;
	
	// Vector of covariate features for each bin
	vector<vector<double> > covar_features;
	
	// Number of clusters to create: numclust
	unsigned int numclust = 100;
	
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
	// return 0;
// 	FILE *testfile_ptr = fopen("test-bin-code/testfile.txt", "w");
// 	for (unsigned int i = 0; i < ann_array.size(); i++) {
// 		fprintf(testfile_ptr, "%s\t%s\t%s\n", ann_array[i][0].c_str(), ann_array[i][1].c_str(), ann_array[i][2].c_str());
// 	}
// 	fclose(testfile_ptr);
// 	return 0;
	
	// Subtract the portions of each annotation that overlap prohibited regions
	for (unsigned int i = 0; i < ann_array.size(); i++) {
	
		// DEBUG
		// printf("Loop 1: %d\n", (int)i);
	
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
	
		// DEBUG
		// printf("Loop 2: %d\n", (int)i);
	
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
		
		// DEBUG
		// printf("Loop 3\n");
		
		ann_array.erase(ann_array.end());
	}
	
	// DEBUG
	// return 0;
	
	// DEBUG - check the genome bins
	// printf("%s\t%s\t%s\n", ann_array[0][0].c_str(), ann_array[0][1].c_str(), ann_array[0][2].c_str());
// 	string debug_file = "/net/gerstein/ll426/code/moat-test/test-bin-code/debug.txt";
// 	FILE *debug_ptr = fopen(debug_file.c_str(), "w");
// 	for (unsigned int i = 0; i < ann_array.size(); i++) {
// 		// printf("%d\n", i);
// 		fprintf(debug_ptr, "%s\t%s\t%s\n", ann_array[i][0].c_str(), ann_array[i][1].c_str(), ann_array[i][2].c_str());
// 	}
// 	fclose(debug_ptr);
	// Early termination
	// return 0;
	
	/* Begin building covariate signal profiles of the genome bins */
	// First step is to produce a file of genome bins
	string regions_presig = "regions.bed";
	int regnum = 1;
	FILE *regions_presig_ptr = fopen(regions_presig.c_str(), "w");
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		char regnum_cstr[STRSIZE];
		sprintf(regnum_cstr, "%d", regnum);
		string outstring = ann_array[i][0] + "\t" + ann_array[i][1] + "\t" + ann_array[i][2] + "\t" + "reg" + string(regnum_cstr) + "\n";
		fprintf(regions_presig_ptr, "%s", outstring.c_str());
		regnum++;
	}
	fclose(regions_presig_ptr);
	
	string regions_postsig = "signals.bed";
	string regions_postsig_sorted = "signals_sorted.txt";
	for (unsigned int i = 0; i < covar_files.size(); i++) {
		string covar_file = covar_files[i];
		
		string command = "./bigWigAverageOverBed " + covar_file + " " + regions_presig + " " + regions_postsig;
		system(command.c_str());
		
		// Next command depends on OS
		command = "uname";
		char buf[STRSIZE];
		string os = "";
		FILE* pipe = popen(command.c_str(), "r");
		if (!pipe) throw runtime_error("Could not determine operating system. Exiting.\n");
		try {
			while (!feof(pipe)) {
				if (fgets(buf, 128, pipe) != NULL) {
					os += buf;
				}
			}
		} catch (...) {
			pclose(pipe);
			throw;
		}
		pclose(pipe);
		
		if (os == "Darwin\n") { // OS == Mac OS X
			command = "sed -i .bak 's/^reg//g' " + regions_postsig;
		} else { // Assume Linux, or other OS that is compatible with this command
			command = "sed -i 's/^reg//g' " + regions_postsig;
		}
		system(command.c_str());
		
		command = "sort -n -k 1,1 " + regions_postsig + " > " + regions_postsig_sorted;
		system(command.c_str());
		
		// DEBUG
		// return 0;
		
		// Read the output into memory to add to the covariate features vector
		int line_index = 0;
		char linebuf_cstr[STRSIZE];
		FILE *regions_postsig_ptr = fopen(regions_postsig_sorted.c_str(), "r");
		while (fgets(linebuf_cstr, STRSIZE-1, regions_postsig_ptr) != NULL) {
		
			string linebuf = string(linebuf_cstr);
			int col_index = 0;
			while (col_index < 5) {
				unsigned int pos = linebuf.find_first_of("\t");
				linebuf = linebuf.substr(pos+1);
				col_index++;
			}
	
			// Now linebuf has the value we're looking for. Put it in the covar_features vector.
			double feature;
			sscanf(linebuf.c_str(), "%lf", &feature);
			
			if (covar_features.size() <= (unsigned int)line_index) {
				vector<double> temp;
				temp.push_back(feature);
				covar_features.push_back(temp);
			} else {
				covar_features[line_index].push_back(feature);
			}
			line_index++;
		}
		if (!(feof(regions_postsig_ptr)) && ferror(regions_postsig_ptr)) { // This is an error
			char preamble[STRSIZE];
			sprintf(preamble, "There was an error reading from %s", regions_postsig_sorted.c_str());
			perror(preamble);
			return 1;
		}
		fclose(regions_postsig_ptr);
	}
	
	// DEBUG - check the covariate vector
// 	for (unsigned int i = 0; i < covar_features.size(); i++) {
// 		for (unsigned int j = 0; j < covar_features[i].size(); j++) {
// 			printf("%f\n", covar_features[i][j]);
// 		}
// 	}
// 	return 0;
		
	/* Cluster the covariate matrix rows */
	
	// Z normalize the data
	for (unsigned int i = 0; i < covar_features[0].size(); i++) {
		double mean = 0.0;
		double var = 0.0;
		
		for (unsigned int j = 0; j < covar_features.size(); j++) {
			mean += covar_features[j][i];
		}
		mean = mean/(double)covar_features.size();
		
		for (unsigned int j = 0; j < covar_features.size(); j++) {
			var += pow((covar_features[j][i]-mean), 2.0);
		}
		var = var/(double)covar_features.size();
		double sd = sqrt(var);
		
		for (unsigned int j = 0; j < covar_features.size(); j++) {
			covar_features[j][i] = (covar_features[j][i]-mean)/sd;
		}
	}
	
	// DEBUG - check scaled values
// 	string scaled_file = "/net/gerstein/ll426/code/moat/scaled.txt";
// 	FILE *scaled_ptr = fopen(scaled_file.c_str(), "w");
// 	for (unsigned int i = 0; i < covar_features.size(); i++) {
// 		fprintf(scaled_ptr, "%f\n", covar_features[i][0]);
// 	}
// 	fclose(scaled_ptr);
// 	return 0;
	
	// Instantiate numclust cluster centroids
	// i.e. Pick numclust rows from the data
	unsigned int step_size = covar_features.size()/numclust;
	vector<vector<double> > centroids;
	for (unsigned int i = 0; i < numclust; i++) {
		centroids.push_back(covar_features[i*step_size]);
	}
	
	// Vector of cluster memberships
	vector<unsigned int> member (covar_features.size(),0);
	
	// Map of the empty clusters
	map<unsigned int,int> empty;
	
	// DEBUG
	// return 0;
	
	// Run for a maximum of 10 iterations
	for (int m = 0; m < 10; m++) {
		// printf("%d\n", m); // DEBUG
		// Calculate cluster memberships
		for (unsigned int i = 0; i < covar_features.size(); i++) {
			vector<double> row = covar_features[i];
			unsigned int minclust;
			double dist = DBL_MAX;
			for (unsigned int j = 0; j < numclust; j++) {
			
				if (empty[j]) {
					continue;
				}
			
				vector<double> cand_centroid = centroids[j];
				double cand_dist = euclidean(row, cand_centroid);
				if (cand_dist < dist) {
					minclust = j;
					dist = cand_dist;
				}
			}
			member[i] = minclust;
		}
		
		// DEBUG
		// return 0;
	
		// Calculate new centroid locations
		for (unsigned int i = 0; i < numclust; i++) {
		
			if (empty[i]) {
				continue;
			}
		
			vector<vector<double> > rows;
			for (unsigned int j = 0; j < covar_features.size(); j++) {
				if (member[j] == i) {
					rows.push_back(covar_features[j]);
				}
			}
			
			if (rows.empty()) {
				empty[i] = 1;
				continue;
			}
			
			// DEBUG
// 			printf("Number of rows: %d\n", (int)rows.size());
// 			printf("Number of columns: %d\n", (int)rows[0].size());
// 			printf("First row/col: %f\n", rows[0][0]);
			// return 0; // DEBUG
			
			// Average calculations
			for (unsigned int j = 0; j < rows[0].size(); j++) {
				double sum = 0.0;
				for (unsigned int l = 0; l < rows.size(); l++) {
					sum += rows[l][j];
				}
				// printf("Sum: %f\n", sum); // DEBUG
				double avg = sum/(double)rows.size();
				// printf("Average: %f\n", avg); // DEBUG
				centroids[i][j] = avg;
			}
		}
	}
	
	// ASSERTION: "centroids" represents the converged cluster centers, and "member"
	// indicates the closest centroid to each vector in "covar_features"
	
	// DEBUG
	// printf("Breakpoint 2\n");
	// return 0;
// 	string centroids_file = "/net/gerstein/ll426/code/moat/centroids.txt";
// 	FILE *centroids_ptr = fopen(centroids_file.c_str(), "w");
// 	for (unsigned int i = 0; i < centroids.size(); i++) {
// 		if (empty[i]) {
// 			continue;
// 		}
// 		for (unsigned int j = 0; j < centroids[i].size(); j++) {
// 			fprintf(centroids_ptr, "%f\n", centroids[i][j]);
// 		}
// 	}
// 	fclose(centroids_ptr);
// 	
	string cluster_file = "/net/gerstein/ll426/code/moat/clusters.txt";
	FILE *cluster_ptr = fopen(cluster_file.c_str(), "w");
	for (unsigned int i = 0; i < member.size(); i++) {
		fprintf(cluster_ptr, "%d\n", member[i]);
	}
	fclose(cluster_ptr);
// 	return 0;
	
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
		
		vector<vector<string> > permuted_set;
		
		for (unsigned int j = 0; j < numclust; j++) {
		
			// DEBUG
			// printf("Permuted set size: %d\n (clust: %d)\n", (int)permuted_set.size(), (int)j);
		
			if (empty[j]) {
				continue;
			}
		
			unsigned int variant_pointer = 0;
		
			vector<vector<string> > cluster_bins;
		
			// Gather up bins in this cluster
			for (unsigned int k = 0; k < ann_array.size(); k++) {
				if (member[k] == j) {
					cluster_bins.push_back(ann_array[k]);
				}
			}
			
			sort(cluster_bins.begin(), cluster_bins.end(), cmpIntervals);
			
			// nt for this cluster
			string concat_nt = "";
			
			// All the input (observed) variants spanning the cluster bins
			// Records locations in "epoch" coordinates
			// Epoch coordinates are 1-based
			// Also includes the trinucleotide context, so we only need to make one pass
			// on the reference genome
			vector<pair<int,string> > obs_var_pos;
			
			// This keeps track of the number of nucleotides in previously observed
			// cluster bins so that we can calculate accurate epoch coordinates
			int epoch_nt = 0;
			
			// Gather up the locations of all confidently mapped trinucleotides (capital letters)
			// Coordinates are for the second letter in the trinucleotide (where the actual mutation is located)
			map<string,vector<int> > local_nt;
			
			vector<char> base; // No treble
			base.push_back('A');
			base.push_back('G');
			base.push_back('C');
			base.push_back('T');
			
			for (int z = 0; z < 4; z++) {
				for (int y = 0; y < 4; y++) {
					for (int x = 0; x > 4; x++) {
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
			
			for (unsigned int l = 0; l < cluster_bins.size(); l++) {
			
				int rand_range_start = atoi(cluster_bins[l][1].c_str());
				int rand_range_end = atoi(cluster_bins[l][2].c_str());
			
				vector<string> rand_range = cluster_bins[l];
			
				pair<unsigned int,unsigned int> range = intersecting_variants(var_array, rand_range, variant_pointer);
				variant_pointer = range.first;
			
				// DEBUG
				// printf("Breakpoint 3\n");
				// printf("%d,%d\n", range.first, range.second);
			
				// int var_subset_count = range.second - range.first + 1;
// 				if (var_subset_count == 0) {
// 					continue;
// 				}
				
				// Populate obs_var_pos
				for (unsigned int m = range.first; m <= range.second; m++) {
					int cur_var_end = atoi(var_array[m][2].c_str());
					int this_epoch = cur_var_end - rand_range_start;
					this_epoch += epoch_nt;
					
// 					stringstream ss;
// 					string cur_nt;
// 					ss << chr_nt[cur_var_end-2];
// 					ss << chr_nt[cur_var_end-1];
// 					ss << chr_nt[cur_var_end];
// 					ss >> cur_nt;
					string placeholder = "";
					
					pair<int,string> variant (this_epoch, placeholder);
					obs_var_pos.push_back(variant);
				}
					
				epoch_nt += (rand_range_end - rand_range_start);
			}
			
			if (obs_var_pos.size() == 0) {
				continue;
			}
			
			// Reset epoch_nt
			epoch_nt = 0;
			
			// Read in reference
			for (unsigned int l = 0; l < cluster_bins.size(); l++) {
				// FASTA import here
				if (last_chr != cluster_bins[l][0]) {
			
					string filename = fasta_dir + "/" + cluster_bins[l][0] + ".fa";
					fasta_ptr = fopen(filename.c_str(), "r");
			
					int first = 1;
					last_chr = cluster_bins[l][0];
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
				
				int rand_range_start = atoi(cluster_bins[l][1].c_str());
				int rand_range_end = atoi(cluster_bins[l][2].c_str());
				
				// Begin indexing
			
				// Save the nt
				concat_nt += chr_nt.substr(rand_range_start, rand_range_end - rand_range_start);
			
				string cur_chr = cluster_bins[l][0];
				for (int k = rand_range_start+1; k <= rand_range_end; k++) { // 1-based index
			
					// Don't read in characters if it will read off either end
					if (k == 1 || k == hg19_coor[cur_chr]) {
						continue;
					}
				
					char nt1 = toupper(chr_nt[k-2]); // 0-based index
					char nt2 = toupper(chr_nt[k-1]); // 0-based index
					char nt3 = toupper(chr_nt[k]); // 0-based index
				
					// Verify there are no invalid characters
					if (nt2 != 'A' && nt2 != 'C' && nt2 != 'G' && nt2 != 'T' && nt2 != 'N') {
						char errstring[STRSIZE];
						sprintf(errstring, "Error: Invalid character detected in FASTA file: %c. Must be one of [AGCTN].\n", nt2);
						printf(errstring);
						return 1;
					}
				
					stringstream ss;
					string cur_nt;
					ss << nt1;
					ss << nt2;
					ss << nt3;
					ss >> cur_nt;
				
					int this_epoch = k - rand_range_start;
					this_epoch += epoch_nt;
					local_nt[cur_nt].push_back(this_epoch);
				}
			
				// End indexing
				epoch_nt += (rand_range_end - rand_range_start);
			}
			
			// Variant processing loop
			for (unsigned int k = 0; k < obs_var_pos.size(); k++) {
			
				// DEBUG
				// printf("%s:%s-%s\n", var_array[k][0].c_str(), var_array[k][1].c_str(), var_array[k][2].c_str());
				// printf("Variant processing loop: iter: %d; clust: %d; perm: %d\n", (int)k, (int)j, i);
			
				// string cur_nt = obs_var_pos[k].second;
				
				char cur_nt1 = toupper(concat_nt[obs_var_pos[k].first-2]);
				char cur_nt2 = toupper(concat_nt[obs_var_pos[k].first-1]); // 0-based index
				char cur_nt3 = toupper(concat_nt[obs_var_pos[k].first]);
				
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
					printf(errstring);
					return 1;
				}
				
				vector<int> pos2;
				for (unsigned int l = 0; l < pos.size(); l++) {
					if (pos[l] != obs_var_pos[k].first) {
						pos2.push_back(pos[l]);
					}
				}
				// Pick new position
				int new_index = rand() % (pos2.size()); // Selection in interval [0,pos2.size()-1]
				
				int new_epoch = pos2[new_index];
				int new_end;
				
				string cluster_chr;
				
				for (unsigned int l = 0; l < cluster_bins.size(); l++) {
					cluster_chr = cluster_bins[l][0];
					int cluster_start = atoi(cluster_bins[l][1].c_str());
					int cluster_end = atoi(cluster_bins[l][2].c_str());
					int cluster_size = (cluster_end - cluster_start);
					
					if (new_epoch > cluster_size) {
						new_epoch -= cluster_size;
					} else {
						// new_end = (new_epoch - cluster_start);
						break;
					}
				}
				
				vector<string> vec;
				vec.push_back(cluster_chr);
				
				char start_cstr[STRSIZE];
				sprintf(start_cstr, "%d", new_end-1); // 0-based
				vec.push_back(string(start_cstr));
				
				char end_cstr[STRSIZE];
				sprintf(end_cstr, "%d", new_end); // 1-based
				vec.push_back(string(end_cstr));
				
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
		
			// last_chr = ann_array[j][0];
		}
		
		// Open new output file
		char perm_num[STRSIZE];
		sprintf(perm_num, "%d", i+1);
		
		string outfile = outdir + "/permutation_" + string(perm_num) + ".txt";
		FILE *outfile_ptr = fopen(outfile.c_str(), "w");
		// vector<int> this_permutation_counts;
		
		for (unsigned int k = 0; k < permuted_set.size(); k++) {
			
			// DEBUG
			// printf("Loop iter: %d; permuted set size: %d; \n", (int)k, (int)permuted_set.size());
		
			fprintf(outfile_ptr, "%s\t%s\t%s\n", permuted_set[k][0].c_str(), permuted_set[k][1].c_str(), permuted_set[k][2].c_str());
		}
		fclose(outfile_ptr);
	}
	
	// DEBUG
	// printf("Breakpoint 3\n");
	
	// Wrap up by removing the temporary files created along the way
// 	string rmcom = "rm " + regions_presig;
// 	system(rmcom.c_str());
// 	rmcom = "rm " + regions_postsig;
// 	system(rmcom.c_str());
// 	rmcom = "rm " + regions_postsig_sorted;
// 	system(rmcom.c_str());
	
	// Verdun
	return 0;
}
