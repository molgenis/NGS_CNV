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

__Usage__
```
python link_bed_to_bam.py \
	-b /path/to/master_bamlist.txt \
	-k /path/to/RDConnect.csv \
	-o /path/to/rd3_data/master_bed_to_bam.txt \
	-r /path/to/rd3_data/rd3_experiment.txt
```


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


## make_master_batches.sh
This small script divides the BAM files per BED file using grep commands. In this, and other documents, I generally refer to each BED file as a batch (as each BED file has a batch of samples) followed by a number. The batch number is based on the order of the BED files in the RDConnect file. So the first file would be named batch1, the second batch2, etc. Afterwards, I manually made the `batch_to_kit.txt` file that links the batch number to the kit name associated with the BED file.


## make_lns.py
This script makes a shell script that links, via softlinks (ln -s), BAM and BAI sample files for a single batch to a specified directory.

__Required parameters__
* sys.argv[1]: Path to a batch sample file listing made by `make_master_batches.sh`
* sys.argv[2]: Path to write the linking script to
* sys.argv[3]: Path to the directory the resulting linking script should place the BAM and BAI files

__Usage__
```
python make_lns.py \
	/path/to/master_batch1.txt \
	/path/to/lnsbatches/lns_batch1.sh \
	/path/to/samples/batch1/
```


## Execute *.sh created by make_lns.py
Simply execute these scripts made by `make_lns.py` to link samples to batch directories.


## determine_ud_sexes.py
Can be used to determine the sexes of the samples that have been labelled UD as claimed sex. To determine the sex of the UD samples, this script requires the ClusterWES coverage files. It collects all X-chromosome coverage values and divides the median by the median of autosomal coverage. If the resulting ratio is smaller than 0.7 the UD sample will be labelled ‘M’ for male. If the ratio is between 0.85 and 1.3 the UD sample will be labelled ‘F’ for female. In all other cases, the UD sample will remain UD. The minimum threshold of 0.85 was chosen as this is also used by the inhouse pipeline.

__Required parameters__
* sys.argv[1]: Table linking sample to sex
* sys.argv[2]: Directory with normalized.coverage.txt files

__Usage__
```
python determine_ud_sexes.py \
	/path/to/bam_to_sex.txt \
	/path/to/solverd/normalized_coverage/
```
