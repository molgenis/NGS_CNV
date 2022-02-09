# CNV calling scripts

__Used abbreviations__
* CAC: CollectAllelicCounts
* CCRS: CallCopyRatioSegments
* CRC: CollectReadCounts
* DRC: DenoiseReadCounts
* MS: ModelSegments
* PDCR: PlotDenoisedCopyRatios
* PMS: PlotModelledSegments
* PON: CreateReadCountPanelOfNormals


## 1: csj3.py
Used to generate the CollectReadCounts jobs for each sample per BED file (batch). Several inputs are required listed below.
Note that for later steps individual scripts were created, described below.

__Required parameters__
* [-no / --]: Number of job scripts to generate
* [-o / --]: Path the job scripts should be written to
* [-g / --]: What type of jobs to generate (in this case 'crc' for CollectReadCounts)
* [-p / --]: Path to the directory that will contain the data of all the GATK CNV calling steps
* [-i / --]: Path to the BAM and BAI sample files of the selected batch
* [-il / --]: Path to the intervallist of the batch (BED file) that has been preprocessed with PreprocessIntervals
* [-jo / --]: Path the CollectReadCounts jobs should write their output to
* [-j / --]: Prefix for the job names
* [--solverd]: Indicate that the jobs are for the Solve-RD project

__Usage__
```
python /home/umcg-mbeukers/scripts/ngs_cnv/script/csj3.py \
	-no 5 \
	-o /path/to/crc_jobs/batch1 \
	-g crc \
	-p /path/to/gatk4_freeze1/ \
	-i /path/to/gatk4_freeze1/samples/batch1 \
	-il /path/to/gatk4_freeze1/ppi/batch1.preprocessed.interval_list \
	-jo /path/to/gatk4_freeze1/crc/batch1/ \
	-j crc_batch1 \
	--solverd
```


## 2: solve_rd_make_pon.py
Used to make CreateReadCountPanelOfNormal jobs for each BED file (batch). Several inputs are required. First, the bam_to_sex.txt file, created with `link_sex_to_bam.py`, is used to be able to know which sample should be placed in with male/female panel of normals. The second is the clusterWES file, containing the merged sample clusters, for the BED file (batch). The third input is the path to the directory containing the read count files produced by the CollectReadCounts step. The last input consists of a label to use for the output (i.e. batch1). 

__Required parameters__
* [-b / --bam-to-sex]: Path to the bam_to_sex.txt file
* [-c / --s4-clusters]: Path to file containing the merged S4 ClusterWES clusters
* [-j / --job-outdir]: Path to the directory the jobs should write their output to
* [-n / --name]: Name of the batch to use as (i.e. batch1)
* [-o / --outdir]: Path other files should be written to
* [-r / --read-counts-dir]: Path to directory containing count files (.hdf5) created by CollectReadCounts jobs
* [-s / --script-outdir]: Path the job scripts should be written to

__Usage__
```
python solve_rd_make_pon.py \
	-b /path/to/bam_to_sex.txt \
	-c /path/to/batch1_clusters_wes.txt \
	-j /path/to/pon/ \
	-n batch1 \
	-o /path/to/samples_to_pon/ \
	-r /path/to/crc/batch1/ \
	-s /path/to/job_scripts/pon_jobs/
```


## 3: get_pon_svd_values.py
Used to collect the SVD values for each panel of normals file in a given directory.

__Required parameters__
* sys.argv[1]: Path to the directory containing panel of normals files (are .hdf5 files)
* sys.argv[2]: Path to output directory to write the SVD values for each panel of normals to

__Usage__
```
python get_pon_svd_values /path/to/pon/batch1/ /path/to/pon_svds/batch1/
```


## 4: plot_svds.R
Used to plot the SVD values, obtained from `get_pon_svd_values.py` for each batch. This scripts contains hard-coded paths, which need to be modified before use, and is ment to be run in R Studio. 


## 5: ud_to_pon.py
Used to determine which UD samples should be denoised with which panel of normals. 

__Required parameters__
* [-c / --s4-cluster]: Combined S4 ClusterWES files
* [-n / --name]: Name of the batch
* [-o / --outdir]: Path to write the files linking UD samples and Panel of normals files to
* [-r / --read-counts-dir]: Path to the directory containing the read counts (output from CRC)
* [-u / --ud-detsex-file]: Path to the file with assigned sexes to UD samples (output from)

__Usage__
```
python ud_to_pon.py \
	-c /path/to/combined_clusters/batch1.txt \
	-n batch1 \
	-o /path/to/pon_ud/ud_to_pon/ \
	-r /path/to/crc/batch1/ \
	-u /path/to/ud_determined_sex.txt
```


## 3: generate_denoise_readcounts.py
Used to generate DenoiseReadCount jobs for each sample for a single BED file (batch). Several inputs are required. First is the number of eigen-samples to denoise the read counts with. This number differs per batch and has previously been determined with `get_pon_svd_values.py` and `plot_svds.R`. Second is the file indicating which sample should be processed with which panel of normals. The other required inputs are the batch label (i.e. batch1) and two paths to output directories, one where the jobscripts should be written to, the other indicating where the DenoiseReadCounts jobs should write their output to.

__Required parameters__
* [-e / --eigen-samples]: The number of eigen-samples that should be used to denoise the read counts
* [-j / --job-out]: Path to the directory the jobs should write their output to
* [-n / --name]: Name of the batch (i.e. batch1)
* [-o / --outdir]: Path the job scripts should be written to
* [-s / --samples-to-pon]: Path to the file indicating which PanelOfNormals should be used for each sample

__Usage__
```
python generate_denoise_readcounts.py \
	-e 6 \
	-j /path/to/drc/batch1 \
	-n batch1 \
	-o /path/to/job_scripts/drc_jobs/batch1 \
	-s /path/to/samples_to_pon/batch1_fsamples.txt
```


## 4: generate_plot_denoised_copy_ratios.py
Used to generate the PlotDenoisedCopyRatio jobs for each samples, for a single BED file (batch). Several inputs are required. First the directory containing both the standardized and denoised subdirectories. Second, a sequence dictionary (.dict file) of the same genome reference used in previous steps such as CollectReadCounts. If such a sequence dictionary doesn't exist it can be created with GATK `CreateSequenceDictionary`. Lastly are the two ouput directories, one indicating where the jobscripts should be written to, the other where the PlotDenoisedCopyRatios jobs should write their output to.

__Required parameters__
* [-d / --drcdir]: 
* [-j / --joboutdir]: 
* [-o / --outdir]: 
* [-s / --sequence-dictionary]: 

__Usage__
```
python generate_plot_denoised_copy_ratios.py \
	-d /path/to/drc/batch1/ \
	-j /path/to/pdcr/batch1/ \
	-o /path/to/job_scripts/pdcr_jobs/batch1/ \
	-s /path/to/genome_reference.dict
```


## 5: generate_collect_allelic_counts.py
Used to generate slurm jobs to perform CollectAllelicCounts for a single BED file (batch). The scripts generates a jobscript for each sample. Three inputs are required. First, is a directory containing BAM + index files (or links to them) for a single BED file. The second input is the genome reference used to map the reads within the BAM files to. The last input is a .vcf.gz file containing common SNPs. See project documentation for more information on the creation of the common SNPs file.

__Required parameters__
* [-i / --indir]: Path to directory with BAM files for a single BED file (batch)
* [-j / --job-outdir]: Path the generated jobs needs to write their output to
* [-o / --outdir]: Path to write the generated jobs_scripts to
* [-p / --prefix]: Prefix to use for the job scripts
* [-r / --genome-reference]: Path to the genome reference fasta file
* [-s / --snp-data]: Path to the common SNPs file

__Usage__
```
python generate_collect_allelic_counts.py \
	-i /path/to/batch1_bams/ \
	-j /path/to/cac/batch1/ \
	-o /path/to/job_scripts/cac_jobs/batch1/ \
	-p cac_batch1 \
	-r /path/to/genome_ref.fa \
	-s /path/to/gnomad_af10_snps.vcf.gz
```


## 6: generate_model_segments.py
Used to generate slurm jobs to perform the ModelSegments step for a single BED file (batch). A job script is generated for each sample. Two inputs are required. First is the path to the directory containing the output from the CollectAllelicCounts step. All .tsv files, containing the allelic counts will be collected and used. Second, is the path to the denoised read counts directory, from which the script will also collect all .tsv files.

__Required parameters__
* [-c / --cac-dir]: Path to the directory containing the allelic counts for the batch
* [-d / --drc-dir]: Path to directory containing the denoised read counts for the batch
* [-o / --outdir]: Path to directory to write the generated scripts to
* [-j / --job-outdir]: Path to directory the jobs should write the output to

__Usage__
```
python generate_model_segments.py \
	-c /path/to/cac/batch1/ \
	-d /path/to/drc/batch1/denoised_counts/ \
	-o /path/to/job_scripts/ms_jobs/ \
	-j /path/to/ms/batch1/
```


## 7: generate_plot_model_segments.py
Used to generate PlotModelledSegments jobs. Several inputs are required, which are mainly path to directories. These are paths to the directories containing the allelic counts, denoised read counts and the modelled segments files. The fourth input is the genome reference sequence dictionary (can be created with `CreateSequenceDictionary`).

__Required parameters__
* [-a / --allelic-dir]: Path to the directory containing allelic count files produced by CollectAllelicCounts
* [-d / --denoised-dir]: Path to the directory containing the denoised read count files produced by DenoiseReadCounts
* [-j / --joboutdir]: Path the PlotModelledSegments jobs should write their output to
* [-m / --modelsegments-dir]: Path to the directory containing the modelled segments files produced by ModelSegments
* [-o / --outdir]: Path to directory to write jobscripts to
* [-s / --sequence-dict]: Path to the genome reference sequence dictionary

__Usage__
```
python generate_plot_model_segments.py \
	- a /path/to/cac/batch1/ \
	- d /path/to/drc/batch1/denoised/ \
	- j /path/to/pms/ \
	- m /path/to/ms/batch1/ \
	- o /path/to/pms_jobs/ \
	- s /path/to/genome_reference.dict
```


## 8: generate_call_copy_ratio_segments.py
Used to generate CallCopyRatioSegments jobs for each sample of a BED file (batch). The path to the folder containing the results of the ModelSegments step should be provided as input. The script will collect the .cr.seg files for each sample which will be used as input for the job.

__Required parameters__
* [-m / --modelsegments-dir]: Path to the directory containing output from the ModelSegments step
* [-j / --joboutdir]: Path the jobs should write their output to
* [-o / --outdir]: Path the job scripts should be written to

__Usage__
```
python generate_call_copy_ratio_segments.py \
	-m /path/to/ms/batch1/ \
	-j /path/to/job_scripts/ccrs_jobs/batch1/ \
	-o /path/to/ccrs/batch1/
```
