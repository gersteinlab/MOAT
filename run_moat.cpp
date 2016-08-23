#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <limits.h>
#include <unistd.h>
#include <vector>
#include <errno.h>

using namespace std;

#define STRSIZE 256

/* This code serves as the front end for users. The user passes all their arguments
 * here, regardless of which version of MOAT they want to run, and this code will
 * set everything up with the appropriate executable on the back end.
 *
 * Type checking performed here
 * Arguments can be in any order
 * [] indicates a list of valid user options
 * CAPITAL_LETTERS indicates a user-supplied number or file
 *
 * Synopsis (MOAT-a): run_moat --algo=a --parallel=[y/n] -n=NUM_PERMUTATIONS 
 * --dmin=MIN_DIST_FOR_RANDOM_BINS --dmax=MAX_DIST_FOR_RANDOM_BINS
 * --blacklist-file=BLACKLIST_FILE --vfile=VARIANT_FILE --afile=ANNOTATION_FILE
 * --out=OUTPUT_FILE
 *
 * Synopsis (MOAT-v): run_moat --algo=v --parallel=[y/n] -n=NUM_PERMUTATIONS 
 * --width=WG_BIN_WIDTH --min_width=MIN_WG_BIN_WIDTH --fasta=WG_FASTA_DIR
 * --blacklist-file=BLACKLIST_FILE --vfile=VARIANT_FILE --out=OUTPUT_DIRECTORY 
 * --ncpu=NUMBER_OF_PARALLEL_CPU_CORES
 *
 * Synopsis (MOATsim): run_moat --algo=s --parallel=[y/n] -n=NUM_PERMUTATIONS 
 * --width=WG_BIN_WIDTH --min_width=MIN_WG_BIN_WIDTH --fasta=WG_FASTA_DIR 
 * --blacklist-file=BLACKLIST_FILE --vfile=VARIANT_FILE --out=OUTPUT_DIRECTORY
 * --ncpu=NUMBER_OF_PARALLEL_CPU_CORES --covar_file=COVARIATE_FILE_1 
 * [--covar_file=COVARIATE_FILE_2 ...]
 */
int main (int argc, char* argv[]) {

	/* User-supplied arguments */
	
	// MOAT algorithm to use
	char algo;
	
	// Use parallel version?
	char parallel;
	
	// Number of permuted variant datasets to create
	int num_permutations = -1;
	
	/* MOAT-a specific arguments */
	
	// The minimum distance between element and random bin
	int dmin = INT_MIN;
	
	// The maximum distance between element and random bin
	int dmax = INT_MIN;
	
	/* MOAT-v specific arguments */
	
	// Width of the bins in which we permute variants
	int width = INT_MIN;
	
	// Minimum width allowed of the variant permutation bins
	int min_width = INT_MIN;
	
	// Directory with wg FASTA files
	string fasta_dir;
	
	// File with prohibited coordinates
	// Expected format: tab(chr, start, end, ...)
	string prohibited_file;
	
	// File with single nucleotide variants
	// Expected format: tab(chr, start, end, ...)
	string vfile;
	
	// File with annotations to study for mutation burden
	// Expected format: tab(chr, start, end, name, ...)
	string afile;
	
	// Output file, or directory with the output files
	// Format: tab(chr, start, end)
	string out;
	
	// Number of CPU cores to use in parallel for MOAT-v
	// If this option is missing, set it to the max number of available CPU cores detected
	// Alternatively, one can use --ncpu=MAX to get the max number of available CPU cores
	int ncpu = INT_MIN;
	
	// Covariate files to use for clustering whole genome bins in MOATsim
	// Must have at least one of these
	vector<string> covar_files;
	
	if (argc >= 10) {
		printf("Usage: run_moat --algo=[a/v/s] --parallel=[y/n] -n=NUM_PERMUTATIONS ");
 		printf("[--dmin=MIN_DIST_FOR_RANDOM_BINS] [--dmax=MAX_DIST_FOR_RANDOM_BINS] ");
 		printf("[--width=WG_BIN_WIDTH] [--min_width=MIN_WG_BIN_WIDTH] [--fasta=WG_FASTA_DIR] ");
 		printf("--blacklist_file=BLACKLIST_FILE --vfile=VARIANT_FILE [--afile=ANNOTATION_FILE] ");
 		printf("--out=[OUTPUT_FILE (if using MOAT-a)/OUTPUT_DIRECTORY (if using MOAT-v)] ");
 		printf("[--ncpu=NUMBER_OF_PARALLEL_CPU_CORES] Exiting.\n");
		return 1;
	}
	
	// Process the input
	// Split into name and value, and store
	for (int i = 1; i < argc; i++) {
		string arg = string(argv[i]);
		size_t pos_equals = arg.find('=');
		if (pos_equals == string::npos) {
			printf("Error: Argument missing value: %s. Exiting.\n", arg.c_str());
			return 1;
		}
		string name = arg.substr(0, pos_equals);
		string value = arg.substr(pos_equals+1);
		
		if (name == "--algo") {
			algo = value[0];
		} else if (name == "--parallel") {
			parallel = value[0];
		} else if (name == "-n") {
			num_permutations = atoi(value.c_str());
		} else if (name == "--dmin") {
			dmin = atoi(value.c_str());
		} else if (name == "--dmax") {
			dmax = atoi(value.c_str());
		} else if (name == "--width") {
			width = atoi(value.c_str());
		} else if (name == "--min_width") {
			min_width = atoi(value.c_str());
		} else if (name == "--fasta") {
			fasta_dir = value;
		} else if (name == "--blacklist_file") {
			prohibited_file = value;
		} else if (name == "--vfile") {
			vfile = value;
		} else if (name == "--afile") {
			afile = value;
		} else if (name == "--out") {
			out = value;
		} else if (name == "--ncpu") {
			if (value == "MAX") {
				ncpu = sysconf( _SC_NPROCESSORS_ONLN );
			} else {
				ncpu = atoi(value.c_str());
			}
		} else if (name == "--covar_file") {
			covar_files.push_back(value);
		} else { // User put in an invalid argument
			printf("Error: Invalid argument name: %s. Exiting.\n", name.c_str());
			return 1;
		}
	}
		
	if (algo == 'a') { // MOAT-a
	
		// Check parallel compatibility
		if (parallel == 'y') {
			string parallel_program = "./moat_a_gpu";
			struct stat buf;
			if (stat(parallel_program.c_str(), &buf)) { // Report the error and exit
				printf("Error: No parallel MOAT-a available. Exiting.\n");
				return 1;
			}
		} else if (parallel != 'n') { // Invalid value for parallel flag
			printf("Error: Invalid value for \"--parallel\": %c. Exiting.\n", parallel);
			return 1;
		}
			
		// Verify all relevant variables are defined
		if (num_permutations == INT_MIN) {
			printf("Error: n is not defined. Exiting.\n");
			return 1;
		} else if (num_permutations <= 0) { // Out of range
			printf("Error: n must be at least 1. Exiting.\n");
			return 1;
		}
		if (dmin == INT_MIN) {
			printf("Error: dmin is not defined. Exiting.\n");
			return 1;
		} else if (dmin < 0) { // Out of range
			printf("Error: dmin cannot be negative. Exiting.\n");
			return 1;
		}
		if (dmax == INT_MIN) {
			printf("Error: dmax is not defined. Exiting.\n");
			return 1;
		} else if (dmax < 0) { // Out of range
			printf("Error: dmax cannot be negative. Exiting.\n");
			return 1;
		}
		if (prohibited_file.empty()) {
			printf("Error: Blacklist file is not defined. Exiting.\n");
			return 1;
		}
		if (vfile.empty()) {
			printf("Error: Variant file is not defined. Exiting.\n");
			return 1;
		}
		if (afile.empty()) {
			printf("Error: Annotation file is not defined. Exiting.\n");
			return 1;
		}
		if (out.empty()) {
			printf("Error: Output file is not defined. Exiting.\n");
			return 1;
		}
		
		char num_permutations_cstr[STRSIZE];
		sprintf(num_permutations_cstr, "%d", num_permutations);
		
		char dmin_cstr[STRSIZE];
		sprintf(dmin_cstr, "%d", dmin);
		
		char dmax_cstr[STRSIZE];
		sprintf(dmax_cstr, "%d", dmax);
		
		string exe;
		if (parallel == 'y') {
			exe = "./moat_a_gpu";
		} else if (parallel == 'n') {
			exe = "./moat_a_cpu";
		}
		
		// execl(exe.c_str(), num_permutations_cstr, dmin_cstr, dmax_cstr, prohibited_file.c_str(), vfile.c_str(), afile.c_str(), out.c_str(), (char *)0);
		string command = exe + " " + string(num_permutations_cstr) + " " + string(dmin_cstr) + " " + string(dmax_cstr) + " " + prohibited_file + " " + vfile + " " + afile + " " + out;
		return system(command.c_str());
	
	} else if (algo == 'v' || algo == 's') { // MOAT-v/MOATsim
		
		// Check parallel compatibility
		if (parallel == 'y') {
			string parallel_program;
			if (algo == 'v') {
				parallel_program = "./moat_v_parallel";
			} else if (algo == 's') {
				parallel_program = "./moatsim_parallel";
			}
			struct stat buf;
			if (stat(parallel_program.c_str(), &buf)) { // Report the error and exit
				if (algo == 'v') {
					printf("Error: No parallel MOAT-v available. Exiting.\n");
				} else if (algo == 's') {
					printf("Error: No parallel MOATsim available. Exiting.\n");
				}
				return 1;
			}
			
			// Check ncpu
			if (ncpu == INT_MIN) {
				ncpu = sysconf( _SC_NPROCESSORS_ONLN );
			} else if (ncpu < 2) {
				if (algo == 'v') {
					printf("Error: Cannot use less than two CPU cores in parallel MOAT-v. Exiting.\n");
				} else if (algo == 's') {
					printf("Error: Cannot use less than two CPU cores in parallel MOATsim. Exiting.\n");
				}
				return 1;
			}
			
		} else if (parallel != 'n') { // Invalid value for parallel flag
			printf("Error: Invalid value for \"--parallel\": %c. Exiting.\n", parallel);
			return 1;
		}
		
		// Verify all relevant variables are defined
		if (num_permutations == INT_MIN) {
			printf("Error: n is not defined. Exiting.\n");
			return 1;
		} else if (num_permutations <= 0) { // Out of range
			printf("Error: n must be at least 1. Exiting.\n");
			return 1;
		}
		if (width == INT_MIN) {
			printf("Error: width is not defined. Exiting.\n");
			return 1;
		} else if (width <= 0) { // Out of range
			printf("Error: width must be positive. Exiting.\n");
			return 1;
		}
		if (min_width == INT_MIN) {
			printf("Error: min_width is not defined. Exiting.\n");
			return 1;
		} else if (min_width <= 0) { // Out of range
			printf("Error: min_width must be positive. Exiting.\n");
			return 1;
		}
		if (fasta_dir.empty()) {
			printf("WG FASTA directory is not defined. Exiting.\n");
			return 1;
		}
		if (prohibited_file.empty()) {
			printf("Error: Blacklist file is not defined. Exiting.\n");
			return 1;
		}
		if (vfile.empty()) {
			printf("Error: Variant file is not defined. Exiting.\n");
			return 1;
		}
		if (out.empty()) {
			printf("Error: Output directory is not defined. Exiting.\n");
			return 1;
		}
		
		// MOATsim check for required covariate files
		if (algo == 's') {
			if (covar_files.size() < 1) {
				printf("Error: Must have at least one covariate signal file for MOATsim. Exiting.\n");
				return 1;
			}
			for (unsigned int i = 0; i < covar_files.size(); i++) {
				struct stat cbuf;
				if (stat(covar_files[i].c_str(), &cbuf)) { // Report the error and exit
					printf("Error trying to stat %s: %s\n", covar_files[i].c_str(), strerror(errno));
					return 1;
				}
			}
		}
		
		char num_permutations_cstr[STRSIZE];
		sprintf(num_permutations_cstr, "%d", num_permutations);
		
		char width_cstr[STRSIZE];
		sprintf(width_cstr, "%d", width);
		
		char min_width_cstr[STRSIZE];
		sprintf(min_width_cstr, "%d", min_width);
		
		string exe;
		if (parallel == 'y') {
		
			char ncpu_cstr[STRSIZE];
			sprintf(ncpu_cstr, "%d", ncpu);
		
			if (algo == 'v') {
				exe = "mpirun -n " + string(ncpu_cstr) + " ./moat_v_parallel";
			} else if (algo == 's') {
				exe = "mpirun -n " + string(ncpu_cstr) + " ./moatsim_parallel";
			}
		} else if (parallel == 'n') {
			if (algo == 'v') {
				exe = "./moat_v_serial";
			} else if (algo == 's') {
				exe = "./moatsim";
			}
		}
		
		// execl(exe.c_str(), num_permutations_cstr, width_cstr, min_width_cstr, prohibited_file.c_str(), fasta_dir.c_str(), vfile.c_str(), out.c_str(), (char *)0);
		string command = exe + " " + string(num_permutations_cstr) + " " + string(width_cstr) + " " + string(min_width_cstr) + " " + prohibited_file + " " + fasta_dir + " " + vfile + " " + out;
		if (algo == 's') {
			for (unsigned int i = 0; i < covar_files.size(); i++) {
				command += " ";
				command += covar_files[i];
			}
		}
		return system(command.c_str());
		
	} else { // Algo has an invalid value
		printf("Error: Invalid value for \"--algo\": %c. Exiting.\n", algo);
		return 1;
	}
	return 0;
}
