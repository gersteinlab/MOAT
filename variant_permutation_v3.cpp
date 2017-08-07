#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <errno.h>
#include <sys/stat.h>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <stdexcept>

using namespace std;

#define STRSIZE 1000

// In v2, the variants are treated as positions on a number line that is shifted
// everywhere in the genome, with wraparound between chromosomes

// In v3, only the variants within some number of basepairs upstream and downstream
// of each annotation are shuffled with random uniform distribution

// Like the name says, this subroutine compares intervals, like those found in
// var_array and ann_array
// Also want to guarantee that annotations are sorted on name too, so if there's
// a fourth element in the vector, we compare that too if all other elements are
// equal
bool cmpIntervals (vector<string> a, vector<string> b) {
	
	string a_chr = a[0];
	int a_start = atoi(a[1].c_str());
	int a_end = atoi(a[2].c_str());
	string a_name;
	
	if (a.size() > 3) {
		a_name = a[3];
	}
	
	string b_chr = b[0];
	int b_start = atoi(b[1].c_str());
	int b_end = atoi(b[2].c_str());
	string b_name;
	
	if (b.size() > 3) {
		b_name = b[3];
	}
	
	int a_chr_int;
	if (a_chr == "chrX") {
		a_chr_int = 23;
	} else if (a_chr == "chrY") {
		a_chr_int = 24;
	} else if (a_chr == "chrM" || a_chr == "chrMT") {
		a_chr_int = 25;
	} else if (a_chr == "chrNo") { // Sentinel value
		a_chr_int = 26;
	} else {
		// printf("Test e: %s, %s, %s\n", a[0].c_str(), a[1].c_str(), a[2].c_str()); // DEBUG
		// printf("Test e\n"); // DEBUG
		a_chr = a_chr.substr(3);
		// printf("Test f\n"); // DEBUG
		a_chr_int = atoi(a_chr.c_str());
	}
	
	int b_chr_int;
	if (b_chr == "chrX") {
		b_chr_int = 23;
	} else if (b_chr == "chrY") {
		b_chr_int = 24;
	} else if (b_chr == "chrM" || b_chr == "chrMT") {
		b_chr_int = 25;
	} else if (b_chr == "chrNo") { // Sentinel value
		b_chr_int = 26;
	} else {
		// printf("Test g\n"); // DEBUG
		b_chr = b_chr.substr(3);
		// printf("Test h\n"); // DEBUG
		b_chr_int = atoi(b_chr.c_str());
	}
	
	if (a_chr_int != b_chr_int) {
		return (a_chr_int < b_chr_int);
	} else if (a_start != b_start) {
		return (a_start < b_start);
	} else if ((a.size() == 3) || (a.size() > 3 && a_end != b_end)) {
		return (a_end < b_end);
	} else if (a.size() > 3) {
		// If we got here, then chr, start, and end are all equal. If a name is present,
		// compare on that
		return a_name < b_name;
	} else { // The intervals are equal, so return false since a is not less than b
		return false;
	}
}

// This subroutine compares intervals, specifically testing whether interval a
// starts past the end of b
bool isGreaterIntervals (vector<string> a, vector<string> b) {
	
	string a_chr = a[0];
	int a_start = atoi(a[1].c_str());
	// int a_end = atoi(a[2].c_str());
	
	string b_chr = b[0];
	// int b_start = atoi(b[1].c_str());
	int b_end = atoi(b[2].c_str());
	
	int a_chr_int;
	if (a_chr == "chrX") {
		a_chr_int = 24;
	} else if (a_chr == "chrY") {
		a_chr_int = 25;
	} else if (a_chr == "chrM" || a_chr == "chrMT") {
		a_chr_int = 26;
	} else {
		// printf("Test e: %s, %s, %s\n", a[0].c_str(), a[1].c_str(), a[2].c_str()); // DEBUG
		// printf("Test e\n"); // DEBUG
		a_chr = a_chr.substr(3);
		// printf("Test f\n"); // DEBUG
		a_chr_int = atoi(a_chr.c_str());
	}
	
	int b_chr_int;
	if (b_chr == "chrX") {
		b_chr_int = 24;
	} else if (b_chr == "chrY") {
		b_chr_int = 25;
	} else if (b_chr == "chrM" || b_chr == "chrMT") {
		b_chr_int = 26;
	} else {
		// printf("Test g\n"); // DEBUG
		b_chr = b_chr.substr(3);
		// printf("Test h\n"); // DEBUG
		b_chr_int = atoi(b_chr.c_str());
	}
	
	if (a_chr_int != b_chr_int) {
		return (a_chr_int > b_chr_int);
	} else {
		return (a_start > b_end);
	}
}

// Subroutine for calculating intersecting intervals
vector<vector<string> > intersecting_intervals (vector<vector<string> > intervals, vector<string> region) {
	
	// Output vector
	vector<vector<string> > output_intervals;
	
	// Unpack region
	string region_chr = region[0];
	string region_start = region[1];
	string region_end = region[2];
	int region_start_num = atoi(region_start.c_str());
	int region_end_num = atoi(region_end.c_str());
	
	// Main loop
	for (unsigned int i = 0; i < intervals.size(); i++) {
		// Unpack interval
		string val_chr = intervals[i][0];
		string val_start = intervals[i][1];
		string val_end = intervals[i][2];
		int val_start_num = atoi(val_start.c_str());
		int val_end_num = atoi(val_end.c_str());
		
		// Intersection test
		if (region_chr == val_chr && region_start_num <= val_end_num && region_end_num >= val_start_num) {
			int max_start = max(region_start_num, val_start_num);
			int min_end = min(region_end_num, val_end_num);
			
			vector<string> vec;
			vec.push_back(region_chr);
		
			char max_start_cstr[STRSIZE];
			sprintf(max_start_cstr, "%d", max_start);
			vec.push_back(string(max_start_cstr));
		
			char min_end_cstr[STRSIZE];
			sprintf(min_end_cstr, "%d", min_end);
			vec.push_back(string(min_end_cstr));
		
			output_intervals.push_back(vec);
		}
			
		// Do we continue?
		if (isGreaterIntervals(intervals[i], region)) {
			break;
		}
	}
	return output_intervals;
}

// Subroutine for calculating intersecting variants
// Differs from intersecting_intervals in that two pointers to the first and last of
// the var_array range are returned, rather than the variants themselves
pair<unsigned int,unsigned int> intersecting_variants (vector<vector<string> > const &intervals, vector<string> region, unsigned int left_pointer) {
	
	// Output vector
	vector<vector<string> > output_intervals;
	
	// Unpack region
	string region_chr = region[0];
	string region_start = region[1];
	string region_end = region[2];
	int region_start_num = atoi(region_start.c_str());
	int region_end_num = atoi(region_end.c_str());
	
	unsigned int right_pointer = left_pointer;
	
	// Main loop
	// vector<vector<string> > variants = *intervals;
	for (unsigned int i = left_pointer; i < intervals.size(); i++) {
		// Unpack interval
		string val_chr = intervals[i][0];
		string val_start = intervals[i][1];
		string val_end = intervals[i][2];
		int val_start_num = atoi(val_start.c_str());
		int val_end_num = atoi(val_end.c_str());
		
		// Intersection test
		if (!(isGreaterIntervals(intervals[i], region))) {
			if (region_chr == val_chr && region_start_num <= val_end_num && region_end_num > val_start_num) {
				right_pointer++;
			} else {
				left_pointer++;
				right_pointer++;
			}
		} else { // Verdun
		
			// DEBUG
// 			printf("val_chr: %s\n", val_chr.c_str());
// 			printf("val_start: %d\n", val_start_num);
// 			printf("val_end: %d\n", val_end_num);
// 			printf("reg_chr: %s\n", region_chr.c_str());
// 			printf("reg_start: %d\n", region_start_num);
// 			printf("reg_end: %d\n", region_end_num);
// 			printf("left: %d\n", (int)left_pointer);
// 			printf("right: %d\n", (int)right_pointer);
		
			break;
		}
	}
	pair<unsigned int,unsigned int> output (left_pointer,right_pointer);
	return output;
}

// Subroutine for calculating intersecting variants (MOAT-v specific)
// Differs from intersecting_intervals in that two pointers to the first and last of
// the var_array range are returned, rather than the variants themselves
pair<unsigned int,unsigned int> intersecting_variantsV (vector<vector<string> > const &intervals, vector<string> region, unsigned int left_pointer) {
	
	// Output vector
	vector<vector<string> > output_intervals;
	
	// Unpack region
	string region_chr = region[0];
	string region_start = region[1];
	string region_end = region[2];
	int region_start_num = atoi(region_start.c_str());
	int region_end_num = atoi(region_end.c_str());
	
	unsigned int right_pointer = left_pointer;
	
	// Main loop
	// vector<vector<string> > variants = *intervals;
	for (unsigned int i = left_pointer; i < intervals.size(); i++) {
		// Unpack interval
		string val_chr = intervals[i][0];
		string val_start = intervals[i][1];
		string val_end = intervals[i][2];
		int val_start_num = atoi(val_start.c_str());
		int val_end_num = atoi(val_end.c_str());
		
		// Intersection test
		if (!(isGreaterIntervals(intervals[i], region))) {
			if (region_chr == val_chr && region_start_num <= val_end_num && region_end_num >= val_start_num) {
				right_pointer++;
			} else {
				left_pointer++;
				right_pointer++;
			}
		} else { // Verdun
		
			// DEBUG
// 			printf("val_chr: %s\n", val_chr.c_str());
// 			printf("val_start: %d\n", val_start_num);
// 			printf("val_end: %d\n", val_end_num);
// 			printf("reg_chr: %s\n", region_chr.c_str());
// 			printf("reg_start: %d\n", region_start_num);
// 			printf("reg_end: %d\n", region_end_num);
// 			printf("left: %d\n", (int)left_pointer);
// 			printf("right: %d\n", (int)right_pointer);
		
			break;
		}
	}
	pair<unsigned int,unsigned int> output (left_pointer,right_pointer);
	return output;
}

// Subtract multiple B regions from a single A region
// Assume the B regions are intersecting intervals with A region, and sorted
vector<vector<string> > subtract_intervals (vector<string> region, vector<vector<string> > intervals) {
	string a_chr = region[0];
	string a_start_str = region[1];
	string a_end_str = region[2];
	
	int a_start = atoi(a_start_str.c_str());
	int a_end = atoi(a_end_str.c_str());
	
	// The difference intervals
	vector<vector<string> > output_intervals;
	
	// Previous end, for tracking purposes
	int previous_end = a_start;
	
	for (unsigned int i = 0; i < intervals.size(); i++) {
		vector<string> b_region = intervals[i];
		string b_chr = b_region[0];
		string b_start_str = b_region[1];
		string b_end_str = b_region[2];
		
		int b_start = atoi(b_start_str.c_str());
		int b_end = atoi(b_end_str.c_str());
		
		// First intersection interval abuts the a_region start
		if (a_start == b_start) {
			previous_end = b_end;
			continue;
		}
		
		string output_chr = a_chr;
		int output_start = previous_end;
		int output_end = b_start;
		
		vector<string> vec;
		vec.push_back(output_chr);
		
		char output_start_cstr[STRSIZE];
		sprintf(output_start_cstr, "%d", output_start);
		vec.push_back(string(output_start_cstr));
		
		char output_end_cstr[STRSIZE];
		sprintf(output_end_cstr, "%d", output_end);
		vec.push_back(string(output_end_cstr));
		
		output_intervals.push_back(vec);
		
		previous_end = b_end;
	}
	
	// Take care of the last one (fencepost)
	if (previous_end != a_end) {
		vector<string> vec;
		vec.push_back(a_chr);
		
		char output_start_cstr[STRSIZE];
		sprintf(output_start_cstr, "%d", previous_end);
		vec.push_back(string(output_start_cstr));
		
		char output_end_cstr[STRSIZE];
		sprintf(output_end_cstr, "%d", a_end);
		vec.push_back(string(output_end_cstr));
		
		output_intervals.push_back(vec);
	}
	
	return output_intervals;
}

// C++ version of the variant permutation code
// The variant array and prohibited regions array are assumed to be sorted when
// passed as arguments
vector<vector<string> > permute_variants (int varcount, vector<string> region) {
	
	// Array for collecting output variants
	vector<vector<string> > out_variants;
	
	// Calculate window length
	int window_length = atoi(region[2].c_str()) - atoi(region[1].c_str());
	
	// Loop through variants, adjusting their positions
	for (int i = 0; i < varcount; i++) {
	
		// DEBUG
		// printf("Permutation variant %d\n", i);
		
		// Pick new position
		int delta = rand() % (window_length); // Selection in interval [0,window_length-1]
		
		vector<string> vec;
		vec.push_back(region[0]);
		
		char new_start_cstr[STRSIZE];
		sprintf(new_start_cstr, "%d", atoi(region[1].c_str()) + delta);
		vec.push_back(string(new_start_cstr));
		
		char new_end_cstr[STRSIZE];
		sprintf(new_end_cstr, "%d", atoi(region[1].c_str()) + delta + 1);
		vec.push_back(string(new_end_cstr));
		
		out_variants.push_back(vec);
		
		// DEBUG
		// printf("Breakpoint 6\n");
	}
	return out_variants;
}

// This method allows the code to run system calls and capture their stdout
string exec (const char* cmd) {
	char buffer[STRSIZE];
	string result = "";
	FILE* pipe = popen(cmd, "r");
	if (!pipe) throw std::runtime_error("popen() failed!");
	try {
  	while (!feof(pipe)) {
  		if (fgets(buffer, STRSIZE, pipe) != NULL) {
  			result += buffer;
  		}
		}
	} catch (...) {
		pclose(pipe);
		throw;
	}
	pclose(pipe);
  return result;
}
