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
#include <stdexcept>
#include "variant_permutation_v3.h"

using namespace std;

#define STRSIZE 1000

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
	
	// Option to specify whether to calculate wg signal scores on the variants
	// 'o': Compute wg signal scores for observed annotations only
	// 'p': Compute wg signal scores for observed and permuted annotations
	// 'n': Do not compute wg signal scores
	char funseq_opt;
	
	// WG signal file to use. Must be bigWig format.
	string signal_file;
	
	if (argc != 9 && argc != 10) {
		fprintf(stderr, "Usage: moat_a_cpu [# of permutations] [d_min] [d_max] [prohibited regions file] [variant file] [annotation file] [output file] [wg signal option (o/p/n)] [wg signal file (optional)]. Exiting.\n");
		return 1;
	} else {
		n = atoi(argv[1]);
		dmin = atoi(argv[2]);
		dmax = atoi(argv[3]);
		prohibited_file = string(argv[4]);
		vfile = string(argv[5]);
		afile = string(argv[6]);
		outfile = string(argv[7]);
		funseq_opt = argv[8][0];
		
		if (funseq_opt != 'o' && funseq_opt != 'p' && funseq_opt != 'n') {
			fprintf(stderr, "Invalid option for wg signal option: \'%c\'. Must be either \'o\' or \'p\' or \'n\'. Exiting.\n", funseq_opt);
			return 1;
		}
		
		if (argc == 10) {
			signal_file = string(argv[9]);
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
		if (stat(prohibited_file.c_str(), &databuf)) { // Report the error and exit
			fprintf(stderr, "Error trying to stat %s: %s\n", prohibited_file.c_str(), strerror(errno));
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
	
	// Merge prohibited regions
	prohibited_regions = merge_intervals(prohibited_regions);
	
	// Remove variants and annotations that intersect the blacklist regions
	vector<vector<string> > var_array_new;
	for (unsigned int i = 0; i < var_array.size(); i++) {
		vector<vector<string> > inter = intersecting_intervals(prohibited_regions, var_array[i]);
		if (inter.size() == 0) {
			var_array_new.push_back(var_array[i]);
		}
	}
	var_array = var_array_new;
	
	vector<vector<string> > ann_array_new;
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		vector<vector<string> > inter = intersecting_intervals(prohibited_regions, ann_array[i]);
		if (inter.size() == 0) {
			ann_array_new.push_back(ann_array[i]);
		}
	}
	ann_array = ann_array_new;
	
	// DEBUG
// 	for (unsigned int k = 0; k < var_array.size(); k++) {
// 		printf("%s, %s, %s\n", var_array[k][0].c_str(), var_array[k][1].c_str(), var_array[k][2].c_str());
// 	}
// 	return 0;

		// wg signal score code
	vector<double> signal_scores;
	// vector<int> signal_overcount (ann_array.size(), 0);
	if (funseq_opt != 'n') {
		
		// Retrieve current working directory for temporary output
		string funseq_outdir = exec("pwd");
		funseq_outdir.erase(funseq_outdir.find_last_not_of(" \n\r\t")+1);
		funseq_outdir += "/tmp";
		
		// Verify that temporary output directory exists, or create it if it doesn't
		string command = "mkdir -p " + funseq_outdir;
		system(command.c_str());
			
		// Produce an input file for bigWigAverageOverBed in the temporary folder
		string avg_infile = funseq_outdir + "/" + "avg_infile.txt";
		string avg_outfile = funseq_outdir + "/" + "avg_outfile.txt";
		int regnum = 1;
		FILE *avg_infile_ptr = fopen(avg_infile.c_str(), "w");
		for (unsigned int i = 0; i < var_array.size(); i++) {
			char regnum_cstr[STRSIZE];
			sprintf(regnum_cstr, "%d", regnum);
			string regnum_str = "reg" + string(regnum_cstr);
			fprintf(avg_infile_ptr, "%s\t%s\t%s\t%s\n", var_array[i][0].c_str(), var_array[i][1].c_str(), var_array[i][2].c_str(), regnum_str.c_str());
			regnum++;
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
		
		// Collect sum of signal scores per annotation
		// vector<vector<string> > signal_output;
		
		// Index to track where we are in the var_array
		unsigned int v_index = 0;
		
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
// 			double signal_score;
// 			sscanf(linebuf.c_str(), "%lf", &signal_score);
			
// 			vector<string> temp;
// 			temp.push_back(var_array[v_index][0]);
// 			temp.push_back(var_array[v_index][1]);
// 			temp.push_back(var_array[v_index][2]);
			var_array[v_index].push_back(linebuf);
			// signal_output.push_back(temp);
			v_index++;
		}
		if (!(feof(avg_outfile_ptr)) && ferror(avg_outfile_ptr)) { // This is an error
			char preamble[STRSIZE];
			sprintf(preamble, "There was an error reading from %s", avg_outfile_sorted.c_str());
			perror(preamble);
			return 1;
		}
		fclose(avg_outfile_ptr);
		
		// Clean up temporary folder
		string rm_com = "rm -rf " + funseq_outdir;
		system(rm_com.c_str());
		
		// Sort
		// sort(signal_output.begin(), signal_output.end(), cmpIntervals);
		
		// Gather up and sum the Funseq values over each annotation
		unsigned int signal_var_pointer = 0;
		for (unsigned int i = 0; i < ann_array.size(); i++) {
			pair<unsigned int,unsigned int> range = intersecting_variants(var_array, ann_array[i], signal_var_pointer);
			signal_var_pointer = range.first;
			double signal_sum = 0.0;
			
			for (unsigned int j = range.first; j < range.second; j++) {
				signal_sum += atof(var_array[j][3].c_str());
			}
			signal_scores.push_back(signal_sum);
		}
	}
	
	// Variables for main loop
	unsigned int var_pointer = 0;
	
	// Collect the output values, will end with same size as ann_array
	vector<double> pvalues;
	
	// Vector of p-values calculated on wg signal score
	vector<double> signal_pvalues;
	
	// Main loop: Iterate through the annotations
	for (unsigned int i = 0; i < ann_array.size(); i++) {
	
		// Unpack the current annotation
		string cur_ann_chr = ann_array[i][0];
		string cur_ann_start = ann_array[i][1];
		string cur_ann_end = ann_array[i][2];
		string cur_ann_name = ann_array[i][3];
		
		// Cast chr, start and end into int
		int cur_ann_chr_num;
		if (cur_ann_chr == "chrX") {
			cur_ann_chr_num = 24;
		} else if (cur_ann_chr == "chrY") {
			cur_ann_chr_num = 25;
		} else if (cur_ann_chr == "chrM") {
			cur_ann_chr_num = 26;
		} else {
			string cur_ann_chr_part = cur_ann_chr.substr(3);
			cur_ann_chr_num = atoi(cur_ann_chr_part.c_str());
		}
		
		int cur_ann_start_num = atoi(cur_ann_start.c_str());
		int cur_ann_end_num = atoi(cur_ann_end.c_str());
		
		// Find the intersecting variants for this annotation k_t
		int target_variants = 0;
		
		// wg signal score variables
		double cur_ann_signal;
		int signal_overcount;
		
		if (funseq_opt == 'p') {
			cur_ann_signal = signal_scores[i];
			signal_overcount = 0;
		}
		
		// Start searching from var_pointer
		// Instantiate cur_var variables
		string cur_var_chr = var_array[var_pointer][0];
		string cur_var_start = var_array[var_pointer][1];
		string cur_var_end = var_array[var_pointer][2];
		
		// Cast chr, start and end into int
		int cur_var_chr_num;
		if (cur_var_chr == "chrX") {
			cur_var_chr_num = 24;
		} else if (cur_var_chr == "chrY") {
			cur_var_chr_num = 25;
		} else if (cur_var_chr == "chrM") {
			cur_var_chr_num = 26;
		} else {
			string cur_var_chr_part = cur_var_chr.substr(3);
			cur_var_chr_num = atoi(cur_var_chr_part.c_str());
		}
		
		int cur_var_start_num = atoi(cur_var_start.c_str());
		int cur_var_end_num = atoi(cur_var_end.c_str());
		
		// Now count the intersecting variants
		// var_pointer points to the "earliest" possible annotation, and vpointer
		// points to the variants up until the last intersecting with the annotation
		unsigned int vpointer = var_pointer;
		
		// DEBUG
		// printf("Breakpoint 2\n");
		// printf("Variant: %s, %s, %s\n", cur_var_chr.c_str(), cur_var_start.c_str(), cur_var_end.c_str());
		// printf("Annotation: %s, %s, %s, %s\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str());
		
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
			if (cur_var_chr == "chrX") {
				cur_var_chr_num = 24;
			} else if (cur_var_chr == "chrY") {
				cur_var_chr_num = 25;
			} else if (cur_var_chr == "chrM") {
				cur_var_chr_num = 26;
			} else {
				string cur_var_chr_part = cur_var_chr.substr(3);
				cur_var_chr_num = atoi(cur_var_chr_part.c_str());
			}
		
			cur_var_start_num = atoi(cur_var_start.c_str());
			cur_var_end_num = atoi(cur_var_end.c_str());
		}
		
		// DEBUG
		// printf("Breakpoint 3\n");
// 		printf("<-- %s, %s, %s, %s -->\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str());
// 		printf("Last variant: %s, %s, %s\n", cur_var_chr.c_str(), cur_var_start.c_str(), cur_var_end.c_str());
// 		printf("Target variants: %d\n", target_variants);
		
		// Pick random bins from surrounding regions
		// We take n/2 from the upstream region, and n/2 from the downstream regions
		// Random number drawn from [0, dmax - dmin - (annotation's length)]
		int range = dmax - dmin - (cur_ann_end_num - cur_ann_start_num + 1);
		vector<vector<string> > upstream_random_bins;
		vector<vector<string> > downstream_random_bins;
		
		// Upstream bin selection
		// Configure where the start of this range is
		string rand_range_chr = cur_ann_chr;
		int rand_range_start = ((cur_ann_start_num + cur_ann_end_num)/2) - dmax;
		
		for (int j = 0; j < n/2; j++) {
			int rand_start = rand() % range;
			
			// Find the bounds of this bin
			string rand_chr = rand_range_chr;
			rand_start += rand_range_start;
			int rand_end = rand_start + (cur_ann_end_num - cur_ann_start_num);
			
			char rand_start_cstr[STRSIZE];
			sprintf(rand_start_cstr,"%d", rand_start);
			string rand_start_str = string(rand_start_cstr);
			char rand_end_cstr[STRSIZE];
			sprintf(rand_end_cstr,"%d",rand_end);
			string rand_end_str = string(rand_end_cstr);
			
			vector<string> vec;
			vec.push_back(rand_chr);
			vec.push_back(rand_start_str);
			vec.push_back(rand_end_str);
			upstream_random_bins.push_back(vec);
		}
		
		// Downstream bin selection
		// Configure where the start of this range is
		rand_range_chr = cur_ann_chr;
		rand_range_start = ((cur_ann_start_num + cur_ann_end_num)/2) + dmin;
		
		for (int j = 0; j < n/2; j++) {
			int rand_start = rand() % range;
			
			// Find the bounds of this bin
			string rand_chr = rand_range_chr;
			rand_start += rand_range_start;
			int rand_end = rand_start + (cur_ann_end_num - cur_ann_start_num);
			
			char rand_start_cstr[STRSIZE];
			sprintf(rand_start_cstr,"%d",rand_start);
			string rand_start_str = string(rand_start_cstr);
			char rand_end_cstr[STRSIZE];
			sprintf(rand_end_cstr,"%d",rand_end);
			string rand_end_str = string(rand_end_cstr);
			
			vector<string> vec;
			vec.push_back(rand_chr);
			vec.push_back(rand_start_str);
			vec.push_back(rand_end_str);
			downstream_random_bins.push_back(vec);
		}
		
		// DEBUG: Print the random bins
// 		for (unsigned int j = 0; j < upstream_random_bins.size(); j++) {
// 			printf("%s, %s, %s\n", upstream_random_bins[j][0].c_str(), upstream_random_bins[j][1].c_str(), upstream_random_bins[j][2].c_str());
// 		}
// 		for (unsigned int j = 0; j < downstream_random_bins.size(); j++) {
// 			printf("%s, %s, %s\n", downstream_random_bins[j][0].c_str(), downstream_random_bins[j][1].c_str(), downstream_random_bins[j][2].c_str());
// 		}
		
		// Find the intersecting variants for the random bins
		// Upstream bins: search backwards from the variant at var_array[var_pointer]
		sort(upstream_random_bins.begin(), upstream_random_bins.end(), cmpIntervals);
		unsigned int vpointer2 = var_pointer;
		
		// A collection of intersecting variants counts from the random bins
		vector<int> varcounts;
		
		// Backwards search!
		unsigned int j = upstream_random_bins.size();
		do {
			j--;
		
			// How many variants intersect this bin?
			int this_variants = 0;
			
			// wg signal score, used only if funseq_opt == 'p'
			double signal_sum = 0;
			
			// Unpack the current annotation
			string upstream_ann_chr = upstream_random_bins[j][0];
			string upstream_ann_start = upstream_random_bins[j][1];
			string upstream_ann_end = upstream_random_bins[j][2];
			
			// Cast to int
			int upstream_ann_chr_num;
			if (upstream_ann_chr == "chrX") {
				upstream_ann_chr_num = 24;
			} else if (upstream_ann_chr == "chrY") {
				upstream_ann_chr_num = 25;
			} else if (upstream_ann_chr == "chrM") {
				upstream_ann_chr_num = 26;
			} else {
				string upstream_ann_chr_part = upstream_ann_chr.substr(3);
				upstream_ann_chr_num = atoi(upstream_ann_chr_part.c_str());
			}
		
			int upstream_ann_start_num = atoi(upstream_ann_start.c_str());
			int upstream_ann_end_num = atoi(upstream_ann_end.c_str());
			
			// Unpack the current variant
			string upstream_var_chr = var_array[vpointer2][0];
			string upstream_var_start = var_array[vpointer2][1];
			string upstream_var_end = var_array[vpointer2][2];
			
			// Cast chr, start and end into int
			int upstream_var_chr_num;
			if (upstream_var_chr == "chrX") {
				upstream_var_chr_num = 24;
			} else if (upstream_var_chr == "chrY") {
				upstream_var_chr_num = 25;
			} else if (upstream_var_chr == "chrM") {
				upstream_var_chr_num = 26;
			} else {
				string upstream_var_chr_part = upstream_var_chr.substr(3);
				upstream_var_chr_num = atoi(upstream_var_chr_part.c_str());
			}
		
			int upstream_var_start_num = atoi(upstream_var_start.c_str());
			int upstream_var_end_num = atoi(upstream_var_end.c_str());
			
			// Now count the intersecting variants
			// vpointer2 points to the "earliest" possible annotation, and vpointer3
			// points to the variants up until the last intersecting with the annotation
			unsigned int vpointer3 = vpointer2;
			
			// While vpointer3 does not go past the current annotation
			while (upstream_var_chr_num > upstream_ann_chr_num || (upstream_var_chr_num == upstream_ann_chr_num && upstream_var_start_num >= upstream_ann_start_num)) {
			
				// If the current variant intersects the current annotation, increment target_variants
				if (upstream_var_chr == upstream_ann_chr && upstream_ann_start_num <= upstream_var_end_num && upstream_var_start_num <= upstream_ann_end_num) {
					this_variants++;
					
					// wg signal score code
					if (funseq_opt == 'p') {
						signal_sum += atof(var_array[vpointer3][3]);
					}
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
			
				upstream_var_chr = var_array[vpointer3][0];
				upstream_var_start = var_array[vpointer3][1];
				upstream_var_end = var_array[vpointer3][2];
			
				// Cast chr, start and end into int
				if (upstream_var_chr == "chrX") {
					upstream_var_chr_num = 24;
				} else if (upstream_var_chr == "chrY") {
					upstream_var_chr_num = 25;
				} else if (upstream_var_chr == "chrM") {
					upstream_var_chr_num = 26;
				} else {
					string upstream_var_chr_part = upstream_var_chr.substr(3);
					upstream_var_chr_num = atoi(upstream_var_chr_part.c_str());
				}
		
				upstream_var_start_num = atoi(upstream_var_start.c_str());
				upstream_var_end_num = atoi(upstream_var_end.c_str());
			}
			
			// this_variants has been settled, save for output
			varcounts.push_back(this_variants);
			
			if (funseq_opt == 'p') {
				if (signal_sum >= cur_ann_signal) {
					signal_overcount++;
				}
			}
		} while (j > 0);
		
		// Downstream bins: a more straight forward search :)
		sort(downstream_random_bins.begin(), downstream_random_bins.end(), cmpIntervals);
		vpointer2 = var_pointer;
		
		for (unsigned int j = 0; j < downstream_random_bins.size(); j++) {
			
			// How many variants intersect this bin?
			int this_variants = 0;
			
			// wg signal score, used only if funseq_opt == 'p'
			double signal_sum = 0;
			
			// Unpack the current annotation
			string downstream_ann_chr = downstream_random_bins[j][0];
			string downstream_ann_start = downstream_random_bins[j][1];
			string downstream_ann_end = downstream_random_bins[j][2];
			
			// Cast to int
			int downstream_ann_chr_num;
			if (downstream_ann_chr == "chrX") {
				downstream_ann_chr_num = 24;
			} else if (downstream_ann_chr == "chrY") {
				downstream_ann_chr_num = 25;
			} else if (downstream_ann_chr == "chrM") {
				downstream_ann_chr_num = 26;
			} else {
				string downstream_ann_chr_part = downstream_ann_chr.substr(3);
				downstream_ann_chr_num = atoi(downstream_ann_chr_part.c_str());
			}
			
			int downstream_ann_start_num = atoi(downstream_ann_start.c_str());
			int downstream_ann_end_num = atoi(downstream_ann_end.c_str());
			
			// Unpack the current variant
			string downstream_var_chr = var_array[vpointer2][0];
			string downstream_var_start = var_array[vpointer2][1];
			string downstream_var_end = var_array[vpointer2][2];
			
			// Cast chr, start and end into int
			int downstream_var_chr_num;
			if (downstream_var_chr == "chrX") {
				downstream_var_chr_num = 24;
			} else if (downstream_var_chr == "chrY") {
				downstream_var_chr_num = 25;
			} else if (downstream_var_chr == "chrM") {
				downstream_var_chr_num = 26;
			} else {
				string downstream_var_chr_part = downstream_var_chr.substr(3);
				downstream_var_chr_num = atoi(downstream_var_chr_part.c_str());
			}
		
			int downstream_var_start_num = atoi(downstream_var_start.c_str());
			int downstream_var_end_num = atoi(downstream_var_end.c_str());
			
			// Now count the intersecting variants
			// vpointer2 points to the "earliest" possible annotation, and vpointer3
			// points to the variants up until the last intersecting with the annotation
			unsigned int vpointer3 = vpointer2;
			
			// While vpointer3 does not go past the current annotation
			while (downstream_var_chr_num < downstream_ann_chr_num || (downstream_var_chr_num == downstream_ann_chr_num && downstream_var_end_num <= downstream_ann_end_num)) {
				
				// If the current variant intersects the current annotation, increment target_variants
				if (downstream_var_chr == downstream_ann_chr && downstream_ann_start_num <= downstream_var_end_num && downstream_var_start_num <= downstream_ann_end_num) {
					this_variants++;
					
					// wg signal score code
					if (funseq_opt == 'p') {
						signal_sum += atof(var_array[vpointer3][3]);
					}
				} else { // Update vpointer2
					if (vpointer3 != var_array.size()-1) {
						vpointer2 = vpointer3 + 1;
					}
				}
				// Now update the cur_var
				if (vpointer3 == var_array.size()-1) {
					break;
				}
				vpointer3++;
				
				downstream_var_chr = var_array[vpointer3][0];
				downstream_var_start = var_array[vpointer3][1];
				downstream_var_end = var_array[vpointer3][2];
			
				// Cast chr, start and end into int
				if (downstream_var_chr == "chrX") {
					downstream_var_chr_num = 24;
				} else if (downstream_var_chr == "chrY") {
					downstream_var_chr_num = 25;
				} else if (downstream_var_chr == "chrM") {
					downstream_var_chr_num = 26;
				} else {
					string downstream_var_chr_part = downstream_var_chr.substr(3);
					downstream_var_chr_num = atoi(downstream_var_chr_part.c_str());
				}
		
				downstream_var_start_num = atoi(downstream_var_start.c_str());
				downstream_var_end_num = atoi(downstream_var_end.c_str());
			}
			
			// this_variants has been settled, save for output
			varcounts.push_back(this_variants);
			
			if (funseq_opt == 'p') {
				if (signal_sum >= cur_ann_signal) {
					signal_overcount++;
				}
			}
		}
		
		// P-value calculation: how many of the random bins have at least as many
		// variants as k_t?
		int overbins = 0;
		for (unsigned int j = 0; j < varcounts.size(); j++) {
			if (varcounts[j] >= target_variants) {
				overbins++;
			}
		}
		
		// DEBUG
		// printf("overbins: %d; n: %d\n", overbins, n);
		
		double fraction = (double)overbins/(double)n;
		pvalues.push_back(fraction);
		
		// DEBUG
		// printf("%d\n", target_variants);
		
		// wg signal p-value calculation
		if (funseq_opt == 'p') {
			double pfrac = (double)signal_overcount/(double)n;
			signal_pvalues.push_back(pfrac);
		}
	}
	
	// Output generation
	FILE *outfile_ptr = fopen(outfile.c_str(), "w");
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		
		// Unpack the annotation
		string cur_ann_chr = ann_array[i][0];
		string cur_ann_start = ann_array[i][1];
		string cur_ann_end = ann_array[i][2];
		string cur_ann_name = ann_array[i][3];
		
		// Print the output line
		if (funseq_opt == 'p') {
			fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%f\t%e\t%f\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str(), pvalues[i], signal_scores[i], signal_pvalues[i]);
		} else if (funseq_opt == 'o') {
			fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%f\t%e\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str(), pvalues[i], signal_scores[i]);
		} else {
			fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%f\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str(), pvalues[i]);
		}
	}
	fclose(outfile_ptr);
	
	// Verdun
	return 0;
}
