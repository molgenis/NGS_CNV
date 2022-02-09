# HC BED file scripts

## generate_sambamba.py
	Script that generates sambamba jobs for a set of BAM files.

## make_hc_scripts.py
	This script will produce modified versions of `make_hcbedfile_part1.sh` and `make_hcbedfile_part2.sh` supplying certain variables with user provided data (such as directory and file paths).
	The following variables will be replaced:
	* $hcdir : Path to directory to use for creating the High Confident BED file.
	* $bedfile : Path to BED file to use as a starting point for creating the High Confident BED file.
	* $slicedbedfile : Path to the sliced BED file (output from `slice_bed_file.py`) that will be further used to make the High Confident BED file.
	* $umcuscriptsdir : Path to the directory containing the UMCU scripts to create the High Confident BED file.
	* $umcuvenvdir : Path to the directory containing the UMNCU Python Virtual Environment (this environment contains several required packages).

## make_hcbedfile_part1.sh
	Template for the first part of creating the High Confident BED file. The first part makes the sliced BED file, the two populations and generates the sambamba jobs.

## make_hcbedfile_part2.sh
	Template for the second part of creating the High Confident BED file. The second part calculates the statistics, merges the populations and eventually creates the High Confident BED file.

## make_population_sample_list.py
	Script to create a BAM list file for the F1, F2, M1 or M2 population

## About the samples list file:
	To make the High Confident BED file, two male and female sample populations of equal sizes are required.
	This samples list file should contain two columns separated by '\t'.
	The first column should be the full path to the BAM files.
	The second column should either be F1, F2, M1 or M2 to denotate which population the sample should be placed in
