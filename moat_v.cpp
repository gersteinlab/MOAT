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
	
	// File with annotations to study for mutation burden
	// Expected format: tab(chr, start, end, name, ...)
	string afile;
	
	// Directory with the output files
	// Format: tab(chr, start, end)
	string outdir;
	
	if (argc != 10) {
		fprintf(stderr, "Usage: moat_v [3mer preservation option (y/n)] [# permuted datasets] [permutation window radius] [min width] [prohibited regions file] [FASTA dir] [variant file] [annotation file] [output directory]. Exiting.\n");
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
		afile = string(argv[8]);
		outdir = string(argv[9]);
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
	
	/* Data structures for the starting data */
	// Variant array, contains variants of the format vector(chr, start, end)
	vector<vector<string> > var_array;
	
	// Annotation array, contains annotations of the format vector(chr, start, end, ann_name)
	// Will be generated from the whole genome coordinates at runtime
	vector<vector<string> > ann_array;
	
	// Prohibited regions array, contains annotations of the format vector(chr, start, end)
	vector<vector<string> > prohibited_regions;
	
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
	// printf("Breakpoint 1\n");
	
	// Sort the arrays
	sort(var_array.begin(), var_array.end(), cmpIntervals);
	sort(ann_array.begin(), ann_array.end(), cmpIntervals);
	sort(prohibited_regions.begin(), prohibited_regions.end(), cmpIntervals);
	
	// Merge annotations and prohibited regions
	ann_array = merge_intervals(ann_array);
	prohibited_regions = merge_intervals(prohibited_regions);
	
	// BEGINNING OF NEW EPOCH CODE
	// Annotation file is relevant again, must import and use in the prohibited
	// region subtraction
	// Removal of hg19 coor code
	
	// Remove each annotation that overlaps prohibited regions
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		vector<vector<string> > int_intervals = intersecting_intervals(prohibited_regions, ann_array[i]);
		
		if (int_intervals.size() > 0) { // Mark for removal
			ann_array[i][0] = "chrNo";
		}
	}
	
	sort(ann_array.begin(), ann_array.end(), cmpIntervals);
	
	// Remove those marked for deletion
	while (ann_array[ann_array.size()-1][0] == "chrNo") {
		ann_array.erase(ann_array.end());
	}
	
	// DEBUG
	// printf("Breakpoint 2\n");
	
	// EPOCH CODE	
	FILE *fasta_ptr = NULL;
	string last_chr = "";
	// int char_pointer;
	string chr_nt;
	int epoch_total = 0;
	string epoch_nt;
	
	// FASTA import and indexing goes here now, before we start permutations
	
	map<string,vector<int> > local_nt;
			
	vector<char> base; // No treble
	base.push_back('A');
	base.push_back('G');
	base.push_back('C');
	base.push_back('T');
	
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
			
// 	if (trimer) {

		unsigned int variant_pointer = 0;
		vector<int> epoch_var;

		for (unsigned int j = 0; j < ann_array.size(); j++) {
			epoch_total += (ann_array[j][2] - ann_array[j][1]);
			if (trimer) {
				if (last_chr != ann_array[j][0]) {
			
					last_chr = ann_array[j][0];
					string chr = last_chr;
					string filename = fasta_dir + "/" + chr + ".fa";
					fasta_ptr = fopen(filename.c_str(), "r");
			
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
						return 1;
					}
				}
			
				epoch_nt += chr_nt.substr(ann_array[j][1], (ann_array[j][2] - ann_array[j][1]));
			}
			
			// Find overlapping variants and convert to epoch coordinates
			// Guarantee the overlapping variants are not going to be less than variant_pointer
			vector<string> ann_range;
			ann_range.push_back(ann_array[j][0]);
			
			char ann_range_start_cstr[STRSIZE];
			sprintf(ann_range_start_cstr, "%d", ann_array[j][1]);
			ann_range.push_back(string(ann_range_start_cstr));
			
			char ann_range_end_cstr[STRSIZE];
			sprintf(ann_range_end_cstr, "%d", ann_array[j][2]);
			ann_range.push_back(string(ann_range_end_cstr));
			
			pair<unsigned int,unsigned int> range = intersecting_variants(var_array, ann_range, variant_pointer);
			variant_pointer = range.first;
			
			for (unsigned int k = range.first; k <= range.second; k++) {
				int this_epoch = epoch_total - (ann_array[j][2] - ann_array[j][1]) + var_array[k][2]; // 1-based
				epoch_var.push_back(this_epoch);
			}
		}
		
		// Index the epoch_nt
		if (trimer) {
			for (int j = 1; j <= (int)epoch_nt.size()-1; j++) {
			// for (int j = 1; j < 20000; j++) { // DEBUG
			
				stringstream ss;
				string cur_nt;
				
				char nt1 = toupper(chr_nt[j-1]);
				char nt2 = toupper(chr_nt[j]);
				char nt3 = toupper(chr_nt[j+1]);
			
				ss << nt1;
				ss << nt2;
				ss << nt3;
				ss >> cur_nt;
				
				// Verify there are no invalid characters
				if (nt2 != 'A' && nt2 != 'C' && nt2 != 'G' && nt2 != 'T' && nt2 != 'N') {
					char errstring[STRSIZE];
					sprintf(errstring, "Error: Invalid character detected in FASTA file: %c. Must be one of [AGCTN].\n", chr_nt[j]);
					fprintf(stderr, errstring);
					return 1;
				}
			
				if (nt1 == 'N' || nt2 == 'N' || nt3 == 'N') {
					continue;
				}
				
					// int this_epoch = epoch_nt + (j+1); // 1-based
				
				local_nt[cur_nt].push_back(j+1); // 1-based
				
				// DEBUG
				// printf("%d\n", this_epoch);
				// printf("%s\n", cur_nt.c_str());
			}
		}
	
	// DEBUG
// 	vector<int> test = local_nt["TCA"];
// 	printf("%d\n", (int)test.size());
// 	for (unsigned int m = 0; m < test.size(); m++) {
// 		printf("%d\n", test[m]);
// 	}
	
	/* Permutate variant locations */
	srand(0);
	
	for (int i = 0; i < num_permutations; i++) {
		
		// Open new output file
		char perm_num[STRSIZE];
		sprintf(perm_num, "%d", i+1);
		
		string outfile = outdir + "/permutation_" + string(perm_num) + ".txt";
		FILE *outfile_ptr = fopen(outfile.c_str(), "w");
	
		// unsigned int variant_pointer = 0;
		
		vector<vector<string> > permuted_set;
			
		// Variant processing loop
		for (unsigned int k = 0; k < epoch_var.size(); k++) {
		
			// DEBUG
			// printf("%s:%s-%s\n", var_array[k][0].c_str(), var_array[k][1].c_str(), var_array[k][2].c_str());
			
			int new_index;
			vector<int> pos2;
			
			// BEGIN 3MER CODE
			if (trimer) {
		
				// vector<string> cur_var = var_array[k];
				char cur_nt1 = toupper(epoch_nt[epoch_var[k]-2]);
				char cur_nt2 = toupper(epoch_nt[epoch_var[k]-1]); // 0-based index
				char cur_nt3 = toupper(epoch_nt[epoch_var[k]]);
			
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
				// pos2 = local_nt[cur_nt];
			
				// If no positions are available, end program with an error and suggest
				// a larger bin size
// 					if (pos.size()-1 == 0) {
// 						char errstring[STRSIZE];
// 						sprintf(errstring, "Error: No valid permutations positions for a variant in bin %s:%s-%s. Consider using a larger bin size.\n",
// 										ann_array[j][0].c_str(), ann_array[j][1].c_str(), ann_array[j][2].c_str());
// 						fprintf(stderr, errstring);
// 						return 1;
// 					}
			
				// vector<int> pos2;
				for (unsigned int l = 0; l < pos.size(); l++) {
					if (pos[l] != epoch_var[k]) {
						pos2.push_back(pos[l]);
					}
				}
				
				if (pos2.size() == 0) {
					continue;
				}
				
				// Pick new position
				new_index = rand() % (pos2.size()); // Selection in interval [0,pos2.size()-1]
				new_index = pos2[new_index];
			} else {
				new_index = rand() % (epoch_total) + 1; // 1-based over whole genome
			}
			
			// END 3MER CODE
			
			vector<string> vec;
			// vec.push_back(var_array[k][0]);
			
			string new_chr;
			
			// EPOCH CODE
			for (unsigned int l = 0; l < ann_array.size(); l++) {
				int ann_size = (ann_array[l][2] - ann_array[l][1]);
				
				if (new_index > ann_size) {
					new_index -= ann_size;
				} else {
					new_chr = ann_array[l][0];
					new_index = ann_array[l][1] + new_index;
					break;
				}
			}
			
			vec.push_back(new_chr);
		
			char start_cstr[STRSIZE];
			sprintf(start_cstr, "%d", new_index-1); // 0-based
			vec.push_back(string(start_cstr));
		
			char end_cstr[STRSIZE];
			sprintf(end_cstr, "%d", new_index); // 1-based
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
		
		for (unsigned int k = 0; k < permuted_set.size(); k++) {
			fprintf(outfile_ptr, "%s\t%s\t%s\n", permuted_set[k][0].c_str(), permuted_set[k][1].c_str(), permuted_set[k][2].c_str());
		}
			// last_chr = ann_array[k][0];
		fclose(outfile_ptr);
	}
	
	// DEBUG
	// printf("Breakpoint 3\n");
	
	// Verdun
	return 0;
}
