# Scripts folder
Contains the scripts used for the CNV calling process and processing of the results. Most folders contain scripts with methods. For example, the `filtering` folder contains several scripts, each script having several methods for a specific filtering such as Conrad filtering.

## Scripts

### classification.py
Classifies GATK4 CNV calls based on the Array CNV calls. There are several classifications: True Positive (TP), False Positive (FP), Array Non-Informative (ANI), WES Non-Informative (WNI) and Array & WES Non-Informatie (AWNI).

__Required parameters__
* [-t / --tool]: The value can be conifer, exomedepth or gatk4.
* [-i / --indir]: Path to the directory with GATK4 CNV calls (CCRS files) to be classified. This parameter should be used to provide GATK4 CNV calls.
* [-if / --infile]: Path to a file containing Conifer/ExomeDepth CNV calls to be classified
* [-a / --arrayfile]: The path to the gold standard file containing the array CNV calls
* [-s / --samples]: Path to the samples table.
* [-e / --exonsfile]: Path to the BED file
* [-p / --probesfile]: Path to the probes file
* [-o / --output]: Path to write the output file to.

__Usage__
```
python scripts/classification.py \
	-t gatk4 \
	-i path/to/gatk4_ccrs_dir \
	-a path/to/goldstandard.txt \
	-e path/to/bedfile.bed \
	-p path/to/array_probes.txt \
	-s path/to/sampletable.txt \
	-o path/to/classification.txt
```


### classify_dualbed.py
Used to classify calls from CNV calling with the normal and the High Confident BED file as shared, overlapping or unique. This produces two new output files, one for the normal and one for the Hiogh-Confident, with the three labels added to the calls. These files I refer to as dualBED.

__Required parameters__
* [-1 / --infile1]: Path to classified CNV calls for the normal BED file
* [-2 / --infile2]: Path to classified CNV calls for the High Confident BED file
* [-o / --outdir]: Path to output directory to write dualBED output files to
* [-p / --percent-overlap]: Minimal percentage overlap required to be considered an overlapping CNV call

__Usage__
```
python classify_dualbed.py \
	-1 /path/to/cnvcalling/classifications.txt \
	-2 /path/to/hc_cnvcalling/hc_classifications.txt \
	-o /path/to/dualbed_classification/ \
	-p 75
```


### comparison.py
Was be used to compare CNV calls from two tools (GATK4 and ExomeDepth for example) that have already been classified with array CNVs.

__Required parameters__
* [-t / --tool]: Type of comparison to perform
* [-1 / --file1]: Path to first list of classified CNV calls
* [-2 / --file2]: Path to second list of CNV calls
* [-l1 / --label1]: Label (name, i.e.: GATK4) for the first list of classified CNVs
* [-l2 / --label2]: Label (name, i.e.: ExomeDepth) for the second list of classified CNVs
* [-o / --outdir]: Path to write output files to
* [-op / --output-prefix]: Prefix for output files

__Usage__
```
python comparison.py \
	-t arraycnvs
	-1 /path/to/cnvcalling/results/gatk4_classification.txt \
	-2 /path/to/cnvcalling/exomedepth/exomedepth_classification.txt \
	-l1 GATK4 \
	-l2 ExomeDepth \
	-o /path/to/cnvcalling/results/comparisons/ \
	-op gatk4_exomedepth
```


### csj2.py
This script was used to generate GATK4 jobs for each CNV calling step as described by the GATK4 tutorial.
Each step has an abbreviation:
* crc = CollectReadCounts
* pon = CreateReadCountPanelOfNormals
* drc = DenoiseReadCounts
* cac = CollectAllelicCounts
* ms = ModelSegments
* ccrs = CallCopyRatioSegments


__Usage CollectReadCounts__
```
python /home/umcg-mbeukers/scripts/ngs_cnv/script/csj2.py \
	-no 28 \
	-o /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/scripts/crc_scripts/28s \
	-g crc \
	-p /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4 \
	-i /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/calling_samples \
	-il /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4 \
	-jo /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/crc/28s \
	-j crc_28s
```

__Usage CreateReadCountPanelOfNormals__
```
python csj2.py \
	-no 1 \
	-o /path/to/cnvcalling/jobs/pon_jobs/
	-g pon
	-p /path/to/cnvcalling/
	-i /path/to/cnvcvalling/reference_counts/
	-jo /path/to/cnvcalling/pon/
	-j reference_pon
```

__Usage DenoiseReadCounts__
```
python /home/umcg-mbeukers/scripts/ngs_cnv/script/csj2.py \
	-no 28 \
	-o /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/scripts/drc_jobs \
	-g drc -p /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4 \
	-i /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/crc/28s \
	-pn /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/pon/18f49m.pon.hdf5 \
	-jo /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/drc \
	-j drc_28s
```

__Usage CollectAllelicCounts__
```
python /home/umcg-mbeukers/scripts/ngs_cnv/script/csj2.py \
	-no 28 \
	-o /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/scripts/cac_jobs \
	-g cac \
	-p /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4 \
	-i /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/calling_samples \
	-r /apps/data/1000G/phase1/human_g1k_v37_phiX.fasta \
	-jo /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/cac \
	-il /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/gnomad_filtering/gnomad_snps.recode.vcf.gz \
	-j cac_28s
```

__Usage ModelSegments__
```
python /home/umcg-mbeukers/scripts/ngs_cnv/script/csj2.py \
	-no 28 \
	-o /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/scripts/ms_jobs \
	-id /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/drc/denoised \
	-ia /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/cac \
	-g ms \
	-p /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4 \
	-i /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4 \
	-jo /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/ms \
	-j ms_28s \
	-t 06:00:00
```

__Usage CallCopyRatioSegments__
```
python /home/umcg-mbeukers/scripts/ngs_cnv/script/csj2.py \
	-no 28 \
	-o /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/scripts/ccrs_scripts \
	-g ccrs \
	-p /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4 \
	-i /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/ms \
	-jo /groups/umcg-gdio/tmp01/umcg-mbeukers/gatk4/ccrs \
	-j ccrs_28s
```


### dualbed.py
Can currently only combine two dualBED files. 

__Required parameters__
* [-t / --tool]: 
* [-1 / --infile1]: Path to the first classified dualBED file
* [-2 / --infile2]: Path to the second classified dualBED file
* [-o / --outfile]: Path to write combined file to

__Usage__
```
python dualbed.py \
	-t combine \
	-1 /path/to/ \
	-2 \
	-o
```

### filtering.py
Scripts that offers several methods of filtering the results. The following filterings can be applied:
* conradfiltering: Filtering CNV calls that have overlapping Conrad CNVs
* nafiltering: Fitering out NAs
* sizefiltering: Filtering out GATK4 CNV calls based on a certain size

#### Conrad filtering
Conrad filtering will filter out GATK4 CNV calls that overlap with a Conrad CNV (which are common CNVs). This filtering can be performed via `python filtering.py -t conradfiltering`.

__Required parameters__
*[-t / --tool]: Type of filtering to perform (in this case ‘conradfiltering’)
*[-i / --infile]: Path to file with classified WES CNV calls
*[-o / --outfile]: Path to write Conrad filtered file to
*[-c / --conradfile]: Path to file containing Conrad CNVs
*[-e / --exonfile]: Path to the BED file to be used (should be the same used for CNV calling)

__Usage__
```
python scripts/filtering.py \
	-t conradfiltering \
	-i path/to/classifications.txt \
	-o path/to/conradfiltered_classifications.txt \
	-c path/to/conradcnvs.txt \
	-e path/to/bedfile.bed
```

#### Size filtering
One of the options is size filtering, which can be done via `python filtering.py -t sizefiltering`. This filters out calls smaller than the provided size.

__Required parameters__
* [-t / --tool]: Type of filtering to perform (in this case ‘sizefiltering’)
* [-i / --infile]: Path to classifications file to filter by call size
* [-o / --outfile]: Path to write filtered output file to
* [-C / --colname]: Name of the column to use for size filtering
* [-S / --cnvsize]: Retain WES CNV calls >= set size

__Usage__
```
python scripts/filtering.py \
	-t sizefiltering \
	-i path/to/classifications.txt \
	-o path/to/classifications_50k.txt \
	-C Conifer_Size \
	-S 50000
```


### totals.py
Can be used to determine the number, and therefore the ratio, of each classification label (True Positive, False Positive, etc) after classification of the CNV calls with the Array calls.

__Required parameters__
* [-t / --tool]: Type of totals to generate (in this case ‘classification’)
* [-i / --infile]: Path to file with classified WES CNV calls
* [-o / --outfile]: Path to write output file to

__Usage__
```
python scripts/totals.py \
	-t classification \
	-i path/to/classifications.txt \
	-o path/to/classification_ratios.txt
```

### utils.py
Can be used to collect some information such as BED regions overlapping with an array CNV, and filtering X&Y chromosomes.

## Folders

### classes/
Contains several python classes used by scripts to store and manipulate data.

### classification/
Contains scripts with methods for the classification of the WES CNV calls with the array calls.

### comparison/
Contains scripts with methods to perform the comparisons.

### filtering/
Contains scripts with methods for the different filtering steps that can be performed with `filtering.py`.

### generate_totals/
Contains scripts for generating total numbers of found array CNVs, False Positives, dualbed ratios, etc.

### gnomad/
Contains scripts to filter the gnomAD file to retain only SNPs with an overal frequency of 10% and to remove indels.

__Usage__
```
python filter_gnomad_af.py \
	-g gnomad.vcf.gz \
	-o gnomad_af10.vcf
```

```
python filter_gnomad_indels.py \
	-g gnomad.vcf.gz \
	-o gnomad_no_indels.vcf
```

### hcbedfile/
Contains two scripts to add genes back to the High Confident BED file. The first is a modified version the UMCUs `slice_bed_file.py`scripts: `alt_slice_bed_file.py` which retains the gene names. The other script `add_genenames_to_hcbedfile.py`adds genenames back to a High Confident BED file.

### parameters/
Contains the scripts that defines the arpgarse parameters for each script as well as methods to check that the required parameters have been set correctly (i.e.: When an input should be a path to a directory it is checked that the provided path actually points to a directory.)

### shared_methods/
Contains a small script with one method most scripts share.

### utils/
Contains scripts for reading input and writing output files as well as several other 
