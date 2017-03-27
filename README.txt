##############################################################################
Mutations Overburdening Annotations Tool (MOAT)
README
v1.0

March 27, 2017

Lucas Lochovsky	and Jing Zhang
Gerstein Lab
Yale University
##############################################################################

Contents:

A) Description
B) Prerequisite Software
C) File List
D) Build Instructions
E) Prerequisite Data Context
F) MOAT-a
G) MOAT-v
H) MOATsim

(A) Description

MOAT (Mutations Overburdening Annotations Tool) is a computational system for identifying significant mutation burdens in genomic elements with an empirical, nonparametric method. Taking a set of variant calls and a set of annotations, both of which may be drawn from anywhere in the human genome, MOAT calculates the observed mutation counts of each annotation, and compares them to the expected mutation counts to detect elevated mutation burdens. The expected mutation count is derived by simulating the expected distribution of background mutations. To produce this expected distribution, MOAT offers two types of permutation algorithm: one that permutes the locations of annotations (MOAT-a), and one that permutes the locations of variants (MOAT-v).

MOAT's annotation permutation algorithm was amenable to parallelization on graphics processing units (GPUs) using Nvidia's CUDA framework due to its high computational intensity and low memory requirements. MOAT's variant permutation algorithm, however, required importing the human reference genome into memory, which made it better suited to parallelization across multiple CPUs with the OpenMPI framework.

Additionally, we provide a variant distribution simulator called MOATsim, which produces permutations of the input variants taking into account trinucleotide identity preservation (similar to MOAT-v), as well as the distribution of whole genome covariates that influence the background mutation rate. Whole genome regions are clustered into groups with similar covariate signal profiles, and intersecting input variants are permuted within all regions in the same cluster, with their new locations preserving the original trinucleotide identity. MOATsim's parallel implementation utilizes the OpenMPI framework.

Furthermore, MOAT offers integration with Funseq2 (funseq2.gersteinlab.org), a framework for evaluating the functional impact of single nucleotide variants. MOAT can compute the Funseq scores of input variants alongside the execution of its permutation algorithms, enabling users to analyze both functional impact and mutation burden simultaneously. MOAT can also compute Funseq scores on the permuted variant datasets generated to simulate the background mutation distribution, resulting in background Funseq scores to compare to the Funseq score of the observed variants.

(B) Prerequisite Software

MOAT is developed for Linux-based operating systems. It is possible that MOAT may function under Mac OS X, or under Windows in an environment like Cygwin, but this has not been tested.

The following software are required to run MOAT. The "version tested" fields indicate the versions of each dependency that have been tested with MOAT. Earlier versions may work, but are unsupported.

1) gcc: The GNU project's C and C++ compiler. Necessary to produce executables from the C++ source code in the MOAT distribution.
	Link: http://gcc.gnu.org/
	Versions tested: 4.4.7, 4.8.2
	
2) make: The GNU make utility, used to automate compilation of C++ source code in MOAT.
	Link: https://www.gnu.org/software/make/
	Version tested: 3.81
	
3) CUDA (For parallel MOAT-a only): Nvidia's Compute Unified Device Architecture framework for parallelizing computational workflows across graphics processing unit (GPU) stream processors. Requires an Nvidia GeForce or Quadro graphics card.
	Link: http://www.nvidia.com/object/cuda_home_new.html
	Version tested: 6.5.12
	
4) OpenMPI (For parallel MOAT-v only): An open source Message Passing Interface implementation that facilitates interprocess communication between parallel processes.
	Link: https://www.open-mpi.org/
	Version tested: 1.10.1

(C) File List

1) bigWigAverageOverBed: Required to compute the covariate signals of each genome bin in MOATsim. Also used to retrieve precomputed Funseq scores. Must be in the MOAT directory.
	*** THE 64-BIT LINUX VERSION OF THIS SCRIPT IS INCLUDED WITH MOAT. FOR OTHER VERSIONS USE THE FOLLOWING LINK.
	Link: http://genome.ucsc.edu/goldenpath/help/bigWig.html (scroll to end of page)

2) makefile: The script that compiles the C++ source code.

3) moat_a.cpp: The C++ source code for the serial version of MOAT-a.

4) moat_a.cu: The source code for the parallel CUDA version of MOAT-a.

5) moat_v.cpp: The C++ source code for the serial version of MOAT-v.

6) moat_v_mpi.cpp: The C++ source code for the parallel OpenMPI version of MOAT-v.

7) moatsim.cpp: The C++ source code for the serial version of MOATsim.

8) moatsim_mpi.cpp: The C++ source code for the parallel OpenMPI version of MOATsim.

9) p_value_emp.cpp: The C++ source code for doing p-value calculations in MOAT-v.

10) README.txt: This file. Contains compilation and usage instructions.

11) run_moat.cpp: The universal frontend for all versions of MOAT. Input checks and data organization are performed here before transferring execution into one of the other executables.

12) variant_permutation_v3.cpp: A collection of helper methods common to all versions of MOAT.

13) variant_permutation_v3.h: Header file corresponding to "variant_permutation_v3.cpp" methods.

(D) Build Instructions

Before MOAT can be used, its C++ source code must be compiled into executable binaries. All the requisite commands have been collected in the "makefile" file in the MOAT directory. To initiate C++ compilation, "cd" into the MOAT directory and run the "make" command.

Command summary:

	cd [MOAT directory]
	make
	
The makefile will compile all the serial source code, and detect any CUDA or OpenMPI installations, compiling the parallel source code if those frameworks are available. The makefile's output will indicate if the parallel versions are available for use.

(E) Prerequisite Data Context

1) hg19 reference genome

MOAT-v's algorithm depends on the human reference genome in order to preserve the trinucleotide context of the input variants when they are permuted. Therefore, MOAT-v requires a local copy of the human genome reference sequence FASTA files to operate. For our analyses, we used the hg19 FASTA files available from the UCSC Genome Browser at:

http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz

2) Funseq whole genome precomputed score file

MOAT-v can compute Funseq scores dynamically using a local installation of Funseq2, or retrieve precomputed scores to greatly accelerate running time (static mode). We provide precomputed Funseq scores on the MOAT website at: []. These scores were generated by running Funseq2 over the entire human genome, saving the highest possible Funseq score at each location. This precomputed score file must be downloaded to the MOAT directory to use the static Funseq analysis of MOAT-v.

(F) MOAT-a

[] = user-defined parameter

Usage: run_moat --algo=a --parallel=[y/n] -n=[number of permutations] --dmin=[minimum distance for random bins] --dmax=[maximum distance for random bins] --blacklist_file=[blacklist file] --vfile=[variant file] --afile=[annotation file] --out=[output file]

MOAT-a is the implementation of the annotation-based permutation algorithm. It iterates through the annotations in the [annotation file], and randomizes [number of permutations (n)] new locations for each annotation in the local genome context, whose boundaries are defined by [minimum distance for random bins (dmin)] and [maximum distance for random bins (dmax)]. Therefore, the local genome context is two intervals: one upstream of the annotation [-dmax, -dmin] and one downstream of the annotation [dmin, dmax]. Variants from the [variant file] are intersected with each of these permuted annotations to produce [n] permuted variant counts. The p-value of each annotation is the fraction of the permuted variant counts that are equal to or greater than the observed variant count (derived from intersecting the [variant file] variants with the annotation). These p-values are written to the [output file].

Additionally, [algo] specifies which permutation algorithm to use: [a] for MOAT-a and [v] for MOAT-v. The [parallel] flag indicates whether to use the CUDA-accelerated version [y] (recommended) or the much slower CPU version [n]. Finally, the [blacklist file] is used to remove regions from consideration that have poor mappability, such as centromeres and telomeres, among others. Annotations that intersect the [blacklist file] will not be analyzed.

(G) MOAT-v

[] = user-defined parameter

Usage: run_moat --algo=v --parallel=[y/n] -n=[number of permutations] --width=[width of whole genome bins] --min_width=[minimum width of whole genome bins] --fasta=[reference genome FASTA file directory] --blacklist_file=[blacklist file] --vfile=[variant file] --out=[output directory] --ncpu=[number of parallel CPU cores to use] --3mer=[y/n]

MOAT-v is the implementation of the variant-based permutation algorithm. It divides the human genome into bins of size [width of whole genome bins], and permutes variants within each bin. Essentially, [width] controls the size of the local genome context, in which we assume the covariates affecting the background mutation rate are more or less constant. The [min_width] parameter influences what MOAT-v will do when it produces bins less than [width], which can happen at the ends of chromosomes, or if a blacklist region is encountered as defined in the [blacklist file]. If the bin's size is higher than [min_width], MOAT-v will proceed with the bin as is. Otherwise, it will attempt to merge the bin with an adjoining neighbor bin to avoid having a bin size below [min_width]. If there are no adjoining neighbors available, the bin will be deleted.

Within each bin, variants from the [variant file] are shuffled to new locations to form [number of permutations (n)] permuted variant datasets. If [3mer] is set to [y]es, these new locations are chosen uniformly over the available trinucleotides that match the trinucleotide identity of the original variant. If [3mer] is set to [n]o, the new locations are chosen uniformly over all bin nucleotides. Trinucleotide identities are derived from the sequence data imported from the [reference genome FASTA file directory]. The [n] permuted variant datasets are produced as [n] files in the [output directory].

Additionally, [algo] specifies which permutation algorithm to use: [a] for MOAT-a and [v] for MOAT-v. The [parallel] flag indicates whether to use the OpenMPI-accelerated version [y] (recommended) or the much slower single CPU version [n]. Finally, the [ncpu] option gives the user control over the number of CPU cores to use in the parallel version. This parameter can be omitted, in which case MOAT-v will automatically use all available CPU cores. This can also be achieved by setting --ncpu to MAX. Alternatively, the user can specify any number between 2 and the maximum number of cores available.

Using run_moat in this way will produce the permuted variant datasets, but to obtain annotation p-values, an additional program must be run: p_value_emp

Usage: p_value_emp [variant file] [annotation file] [prohibited regions file] [permutation variants' directory] [output file]

[variant file] is the same variant file used in run_moat, and the [annotation file] contains the annotations whose mutation burdens you want to evaluate. [prohibited regions file] works the same as [blacklist file], and [permutation variants' directory] corresponds to the [output directory] used in run_moat. The p-values are written to the [output file].

H) MOATsim

[] = user-defined parameter

Usage: run_moat --algo=s --parallel=[y/n] -n=[number of permutations] --width=[width of whole genome bins] --min_width=[minimum width of whole genome bins] --fasta=[reference genome FASTA file directory] --blacklist_file=[blacklist file] --vfile=[variant file] --out=[output directory] --ncpu=[number of parallel CPU cores to use] --3mer=[y/n] --covar_file=[covariate signal file 1] [--covar_file=[covariate_signal_file_2] ...]

MOATsim is a somatic variant simulator that produces [n] permutations of the input [variant file]. Like MOAT-v, it divides the human genome into bins of size [width of whole genome bins]. [width] controls the size of the local genome context, in which we assume the covariates affecting the background mutation rate are more or less constant. The [min_width] parameter influences what MOATsim will do when it produces bins less than [width], which can happen at the ends of chromosomes, or if a blacklist region is encountered as defined in the [blacklist file]. If the bin's size is higher than [min_width], MOATsim will proceed with the bin as is. Otherwise, it will attempt to merge the bin with an adjoining neighbor bin to avoid having a bin size below [min_width]. If there are no adjoining neighbors available, the bin will be deleted.

Bins are clustered using k-means clustering based on the similarity of their covariate signal profiles derived from the [covariate signal files]. At least one such file must be provided when invoking MOATsim. If [3mer] is set to [y]es, the variants from the [variant file] are relocated within their cluster's bins to locations with the same trinucleotide identity as their original sites. If [3mer] is set to [n]o, new locations are chosen uniformly over all nucleotides in the bin cluster. Trinucleotide identities are derived from the sequence data imported from the [reference genome FASTA file directory]. The [n] permuted variant datasets are produced as [n] files in the [output directory].

The [parallel] flag indicates whether to use the OpenMPI-accelerated version [y] (recommended) or the much slower single CPU version [n]. Finally, the [ncpu] option gives the user control over the number of CPU cores to use in the parallel version. This parameter can be omitted, in which case MOATsim will automatically use all available CPU cores. This can also be achieved by setting --ncpu to MAX. Alternatively, the user can specify any number between 2 and the maximum number of cores available.

NOTE: MOATsim relies on the "bigWigAverageOverBed" program to generate the covariate signal profiles for the whole genome bins. The 64-bit Linux version is provided with the MOAT distribution, but if you need another version, other versions are available from:

http://genome.ucsc.edu/goldenpath/help/bigWig.html (scroll to end of page)

Additionally, the use of "bigWigAverageOverBed" will result in the generation of temporary files in the MOAT directory. Therefore, if you attempt to run multiple instances of MOATsim, this will likely cause the data interference between the instances, thereby corrupting the results. Naturally, we recommend you not do this. ;)
