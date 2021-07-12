\scripts
	[make_hcbedfile_part1.sh]: 
	[make_hcbedfile_part2.sh]: 
	[make_hc_scripts.py]: 


\make_hcbedfile_part1.sh
	Template for the first part of creating the High Confident BED file. The first part makes the sliced BED file, the two populations and generates the sambamba jobs.

\make_hcbedfile_part2.sh
	Template fo rthe second part of creating the High Confident BED file. The second part calculates the statistics, merges the populations and eventually creates the High Confident BED file.

\make_hc_scripts.py
	This script will produce modified versions of `make_hcbedfile_part1.sh` and `make_hcbedfile_part2.sh` supplying certain variables with user provided data (such as directory and file paths).
	The following variables will be replaced:
	* $hcdir
	* $bedfile
	* $slicedbedfile

========================================================================================================================================================================================================

About the samples list file:
	To make the High Confident BED file, two male and female sample populations of equal sizes are required.
	This samples list file should contain two columns separated by '\t'.
	The first column should be the full path to the BAM files.
	The second column should either be F1, F2, M1 or M2 to denotate which population the sample should be placed in
