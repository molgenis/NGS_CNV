# Preprocessing scripts

## make_master_bamlist.py
Prior to using this script I made a directory listing of the directory containing the Solve-RD Freeze1 BAM files using the `ls -las` command.
This simple script takes that file and outputs only the full paths to the bam and index files, thus removing any other data outputted by the `ls -las` command. I wrote this output to `master_bamlist.txt`


## link_bed_to_bam.py
Links the BAM files in `master_bamlist.txt` to the appropriate BED files. The resulting output file has two tab separated columns. The first contains the BAM/BAI file. So each sample has two lines, the first with the full path to the BAM file, the second line with the full path to the BAM index file.

__Required parameters__
* [-b / --bamlist]: The `master_bamlist.txt` file
* [-k / --kitlist]: The RDConnect file linking kit name to BED file
* [-o / --outfile]: Path to write `master_bed_to_bam.txt` to
* [-r / --rd3data]: RD3 Experiment table data



## make_master_batches.sh
This small script divides the BAM files per BED file using grep commands. In this, and other documents, I generally refer to each BED file as a batch (as each BED file has a batch of samples) followed by a number. The batch number is based on the order of the BED files in the RDConnect file. So the first file would be named batch1, the second batch2, etc. Afterwards, I manually made the `batch_to_kit.txt` file that links the batch number to the kit name associated with the BED file.


## link_sex_to_bam.py
Links the reported sexes of samples to BAM files using several RD3 data tables. The resulting output file is a three column tab separated table with a header. The three columns are the sample E-number, path to the BAM file and sex of the sample (which is F, M or UD). If the BAM file for a sample was not present, the value will be ‘NA’.

__Required parameters__
* [-b / --bed-to-bam]: Path to the `master_bed_to_bam.txt` file
* [-e / --experiment-data]: Path to the RD3 experiment file
* [-o / --outfile]: Path to write the `sex_to_bam.txt` file to
* [-s / --samples-data]: Path to the RD3 samples file
* [-u / --subjects-data]: Path to the RD3 subjects file

__Usage__
```
python link_sex_to_bam.py \
	-b /path/to/master_bed_to_bam.txt \
	-e /path/to/rd3_data/rd3_experiment.txt \
	-o /path/to/write/sex_to_bam.txt \
	-s /path/to/rd3_data/rd3_samples.txt \
	-u /path/to/rd3_data/rd3_subjects.txt
```


## determine_ud_sexes.py
Can be used to determine the sexes of the samples that have been labelled UD as claimed sex. To determine the sex of the UD samples, this script requires the ClusterWES coverage files. It collects all X-chromosome coverage values and divides the median by the median of autosomal coverage. If the resulting ratio is smaller than 0.7 the UD sample will be labelled ‘M’ for male. If the ratio is between 0.85 and 1.3 the UD sample will be labelled ‘F’ for female. In all other cases, the UD sample will remain UD. The minimum threshold of 0.85 was chosen as this is also used by the inhouse pipeline. This script requires two command line values:
sys.argv[1]: Table linking sample to sex
sys.argv[2]: Directory with normalized.coverage.txt files


## fm_combined_clusters.py
Makes male and female panel of normals for each combined S4 ClusterWES cluster, per batch. First, S4 ClusterWES clusters were combined so each cluster had at least 200 samples. 

__Required parameters__
* [-b / --bam-to-sex]: Path to `sex_to_bam.txt` file created by link_sex_to_bam.py
* [-c / --combined-s4]: Path to combined S4 ClusterWES cluster file
* [-o / --outdir]: Path to write the PoN job scripts to
* [-p / --pon-out]: Path the job should write the PanelOfNormals file to
* [-r / --read-counts-dir]: Path to directory with read count files created by CollectReadCounts

__Usage__
```
python fm_combined_clusters.py
	- 
	- 
```
