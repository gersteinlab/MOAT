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
#include <dirent.h>
#include <stdexcept>
#include "variant_permutation_v3.h"

using namespace std;

#define STRSIZE 256
#define BIGSTRSIZE 10240

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

// This code comes after the empirical mutation permutation code that randomizes
// variant locations. The output from that step is used to calculate the p-values
// of an annotation set given an observed variant set and the permutation variant
// sets

// Variant sets are assumed to be prefiltered to remove those that fall into
// prohibited regions, but the annotation set is assumed not to be. Annotations
// are filtered to remove those in prohibited regions.

int main (int argc, char* argv[]) {

	// File with single nucleotide variants
	// Expected format: tab(chr, start, end, ...)
	string vfile;
	
	// File with annotations to study for mutation burden
	// Expected format: tab(chr, start, end, name, ...)
	string afile;
	
	// File with prohibited coordinates
	// Expected format: tab(chr, start, end, ...)
	string prohibited_file;
	
	// Directory with the permutation variant files. Will parse every file.
	// Expected format: tab(chr, start, end)
	string permutation_dir;
	
	// File with the output
	// Format: tab(chr, start, end, name, p-value)
	string outfile;
	
	// WG signal mode to use
	// 'o': Evaluate wg signal on (o)bserved variants only
	// 'p': Evaluate wg signal on permuted variants as well, output the p-value significance
	// 'n': Do not use wg signal computation
	char funseq_opt;
	
	// WG signal file to use. Must be bigWig format.
	string signal_file;
	
	if (argc != 7 || argc != 8) {
		fprintf(stderr, "Usage: p_value_emp [variant file] [annotation file] [prohibited regions file] [permutation variants' directory] [output file] [wg signal option (o/p/n)] [signal file]. Exiting.\n");
		return 1;
	} else {
		vfile = string(argv[1]);
		afile = string(argv[2]);
		prohibited_file = string(argv[3]);
		permutation_dir = string(argv[4]);
		outfile = string(argv[5]);
		funseq_opt = argv[6][0];
		
		if (argc == 8) {
			signal_file = string(argv[7]);
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
	
	// Verify that "funseq_opt" is valid
	if (funseq_opt != 'o' && funseq_opt != 'p' && funseq_opt != 'n') {
		fprintf(stderr, "Error: Funseq option was set to \'%c\', which is invalid. ", funseq_opt);
		fprintf(stderr, "Must be either \'o\' or \'p\' or \'n\'. Exiting.\n");
		return 1;
	}
	
	// Check bigWigAverageOverBed and Funseq data file in static mode
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
	
	// Collect the output values, will end with same size as ann_array
	vector<double> pvalues;
	
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
	
	vector<int> observed_counts = mutation_counts(var_array, ann_array);
	
	// Vector keeping track of number of permutations where the varcount of an
	// annotation exceeds the observed varcount
	vector<int> num_overcount (ann_array.size(),0);
	
	// Number of permutation files found
	int num_permutations = 0;
	
	DIR *perm_dir;
	struct dirent *ent;
	if ((perm_dir = opendir(permutation_dir.c_str())) != NULL) {
	
		while ((ent = readdir(perm_dir)) != NULL) {
			if (ent->d_type == DT_REG) { // Proceed to open this file and copy the variants into memory
			
				vector<vector<string> > perm_var_array;
				
				// Create the fully qualified filename
				string full_filename;
				if (permutation_dir.find_last_of("/") == permutation_dir.size()-1) { // Ends with a "/"
					full_filename = permutation_dir + string(ent->d_name);
				} else { // Then add the "/"
					full_filename = permutation_dir + "/" + string(ent->d_name);
				}
				
				// Save the first 3 columns, ignore the rest if there are any
				FILE *perm_file = fopen(full_filename.c_str(), "r");
				while (fgets(linebuf, STRSIZE, perm_file) != NULL) {
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
		
					perm_var_array.push_back(vec);
				}
				// Check feof of perm_file
				if (feof(perm_file)) { // We're good
					fclose(perm_file);
				} else { // It's an error
					char errstring[STRSIZE];
					sprintf(errstring, "Error reading from %s", ent->d_name);
					perror(errstring);
					return 1;
				}
				
				vector<int> perm_counts = mutation_counts(perm_var_array, ann_array);
				
				for (unsigned int i = 0; i < ann_array.size(); i++) {
					if (perm_counts[i] >= observed_counts[i]) {
						num_overcount[i]++;
					}
				}
				num_permutations++;
			}
		}
		
	} else {
		// Could not open directory
		char errstring[STRSIZE];
		sprintf(errstring, "Error opening directory %s", permutation_dir.c_str());
		perror(errstring);
		return 1;
	}
	
	// Now do final p-value calculation
	for (unsigned int i = 0; i < ann_array.size(); i++) {
		double fraction = (double)num_overcount[i]/(double)num_permutations;
		pvalues.push_back(fraction);
	}
	
	// Also do wg signal calculation, if requested
	// Legacy code assumes Funseq, hence the variable names
	vector<double> funseq_scores;
	vector<int> fs_overcount (ann_array.size(), 0);
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
		
		// Clean up Funseq temporary folder
		string rm_com = "rm -rf " + funseq_outdir;
		system(rm_com.c_str());
		
		// Additional steps for using permutation Funseq scores
		if (funseq_opt == 'p') {
			for (int i = 1; i <= num_permutations; i++) {
			
				// DEBUG
				// printf("Permutation: %d\n", i);
			
				// Collect sum of Funseq scores per annotation
				vector<double> perm_funseq_scores;
				vector<vector<string> > perm_funseq_output;
				
				char i_str[STRSIZE];
				sprintf(i_str, "%d", i);
			
				// Read in "Output.bed"
				// int first = 1;
				char linebuf3[BIGSTRSIZE];
				string pfunseq_output_file = permutation_dir + "/permutation_" + string(i_str) + ".txt";
				FILE *pffile_ptr = fopen(pfunseq_output_file.c_str(), "r");
				while (fgets(linebuf3, BIGSTRSIZE, pffile_ptr) != NULL) {
				
// 					if (first) {
// 						first = 0;
// 						continue;
// 					}
		
					string line = string(linebuf3);
				
					vector<string> vec;
					for (int i = 0; i < 6; i++) {
						size_t ws_index = line.find_first_of("\t\n");
						string in = line.substr(0, ws_index);
						vec.push_back(in);
						line = line.substr(ws_index+1);
					}
			
					// If this is not a standard chromosome, then remove this row
					if (chr2int(vec[0]) == 0) {
						continue;
					}
			
					perm_funseq_output.push_back(vec);
				}
				// Check feof of vfile
				if (feof(pffile_ptr)) { // We're good
					fclose(pffile_ptr);
				} else { // It's an error
					char errstring[STRSIZE];
					sprintf(errstring, "Error reading from %s", pfunseq_output_file.c_str());
					perror(errstring);
					return 1;
				}
				
				// DEBUG
				// printf("Input step done\n");
			
				// Sort
				sort(perm_funseq_output.begin(), perm_funseq_output.end(), cmpIntervals);
			
				// Gather up and sum the Funseq values over each annotation
				unsigned int pfunseq_var_pointer = 0;
				for (unsigned int j = 0; j < ann_array.size(); j++) {
					pair<unsigned int,unsigned int> range = intersecting_variants(perm_funseq_output, ann_array[j], pfunseq_var_pointer);
					pfunseq_var_pointer = range.first;
					double funseq_sum = 0.0;
				
					for (unsigned int k = range.first; k < range.second; k++) {
						funseq_sum += atof(perm_funseq_output[k][5].c_str());
					}
				
					perm_funseq_scores.push_back(funseq_sum);
				}
				
				// DEBUG
				// printf("Funseq score sum done\n");
				
				// Now update fs_overcount
				for (unsigned int j = 0; j < ann_array.size(); j++) {
				
					// DEBUG
					// printf("Permutation %d, annotation %d: %e, %e\n", i, j, perm_funseq_scores[j], funseq_scores[j]);
				
					if (perm_funseq_scores[j] >= funseq_scores[j]) {
						fs_overcount[j]++;
					}
				}
			}
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
			
			// Now do final p-value calculation
			double fraction = (double)fs_overcount[i]/(double)num_permutations;
		
			fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%f\t%e\t%f\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str(), pvalues[i], funseq_scores[i], fraction);
		} else if (funseq_opt == 'o') {
			fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%f\t%e\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str(), pvalues[i], funseq_scores[i]);
		} else {
			fprintf(outfile_ptr, "%s\t%s\t%s\t%s\t%f\n", cur_ann_chr.c_str(), cur_ann_start.c_str(), cur_ann_end.c_str(), cur_ann_name.c_str(), pvalues[i]);
		}
	}
	fclose(outfile_ptr);
	
	// Verdun
	return 0;
}
