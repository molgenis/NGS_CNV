# Result processing: Filtering scripts

## ccrscall.py
Contains the CcrsCall class that the different filtering scripts use to store GATK4 CNV call data and other data such as overlapping ERN genes.


## combine_calls_per_ern.py
Collects and writes all GATK4 CNV calls per ERN to file. Three inputs are required. The first is a directory containing the combined CCRS files (these are the files in which all separate .called.seg files for samples have been combined). The second input is the path to the directory containing the ERN gene lists. The last input is the path to the `samples_to_ern.txt` file that links samples to ERN. The output is a set of files, one for each ERN containing calls of samples belonging to ERNs. The output files are similar to the combined CCRS files.

__Required parameters__
* [-c / --ccrs-dir]: Path to directory containing the combined CCRS files for all batches
* [-o / --outdir]: Path to directory to write the calls per ERN to
* [-e / --ern-dir]: Path to directory containing the ERN gene lists
* [-s / --samples-to-ern]: Path to `samples_to_ern.txt` file

__Usage__
```
python combine_calls_per_ern.py \
	-c /path/to/combined_ccrs_files/ \
	-o /path/to/calls_per_ern/ \
	-e /path/to/ern_gene_lists/ \
	-s /path/to/samples_to_ern.txt
```


## combine_nonud_ud.py
As samples that had 'UD' as the noted sex were processed in separate folders, this script was used to combine the CNV calls of non-UD (male/female) and UD samples into one single file, one for each BED file (batch). The first input is the directory containing the combined CCRS files (output from `combine_ccrs.py`) for all non-UD samples. The second input is the directory containing the combined CCRS files (output from `combine_ccrs.py`) for all UD samples. The output directory will contain the non-UD and UD CNV calls per BED file (batch).

__Required parameters__
* [-1 / --nonud-dir]: 
* [-2 / --ud-dir]: 
* [-o / --output-dir]: 

__Usage__
```
python comnbine_nonud_ud.py \
	-1 /path/to/combined_ccrs/batch1.called.seg \
	-2 /path/to/ud/combined_ccrs/batch1.called.seg \
	-o /path/to/nonud_ud_combined/
```


## compare_results.py
Used to determine which GATK4 CCRS calls overlap with CNV calls from other tools (ClinCNV, Conifer, ExomeDepth, VarGenius). The first required input is the path to a file containing all GATK4 CNV calls for a specific ERN. This is an output file created by `combine_calls_per_ern.py`. The second input is the path to the file containg CNV calls of one of the other tools (ClinCNV, Conifer, ExomeDepth, VarGenius). These are available on the cluster. The third input is the name of the Solve-RD tool the Solve-RD calls belong to. The value should be `clincnv`, `conifer`, `exomedepth`, or `vargenius`. The fourth input is the ERN gene list, which should be the genelist of the same ERN as the comparison being performed.

__Required parameters__
* [-g / --gatk4-calls]: Path to GATK4 CNV calls for a single ERN
* [-c / --calls]: Solve-RD (ClinCNV, Conifer, ExomeDepth or VarGenius) calls for the same ERN 
* [-ct / --callstool]: Name of the tool that made the calls (can be: clincnv, conifer, exomedepth, vargenius) supplied with the -c/--calls parameter.
* [-o / --outdir]: Path to directory to write output files to
* [-op / --output-prefix]: Output prefix to use
* [-e / --ern-gene-list]: Path to the sepcific ERN genelist

__Usage__
```
python compare_results.py
	-g /path/to/gatk4_rnd.txt \
	-c /path/to/exomedepth_rnd.txt \
	-ct exomedepth \
	-o /path/to/outdir/ \
	-op gatk_exomedepth \
	-e /path/to/ern_genelists/ern_rnd.tsv
```


## conrad_filtering.py
Filters a combined CCRS calls file using Conrad CNV data. For each CCRS call it is determined whether there are overlapping Conrad CNVs and if those Conrad CNVs all have the same call type (duplication or deletion). If this is the case, all gene names overlapping with the CCRS call must also be in the list of gene names associated with the overlapping Conrad CNVs. If that is the case, the call is removed. In all other cases the call is kept.

__Required parameters__
* [-i / --infile]: Path to combined CCRS calls file to filter
* [-o / --outdir]: Path to directory to write output files to
* [-c / --conrad-cnvs]: Path to text file containing the Conrad CNVs
* [-e / --exons-file]: Path to Conrad exons file for autosomal chromosomes
* [-x / --xchrom-exons-file]: Path to Conrad exons file for X-chromosome
* [-p / --out-prefix]: Prefix to use for output files

__Usage__
```
python conrad_filtering.py \
	-i /path/to/combined_ccrs/batch2.called.igv.seg \
	-o /path/to/conrad_filtered/output \
	-c /path/to/conrad_cnvs.txt \
	-e /path/to/conrad_exonshg19.txt \
	-x /path/to/conrad_exonshg19X.txt \
	-p batch2
```


## ern_filtering.py
Can be used to perform the Combined ERN and Singular ERN filtering steps. Calls that overlap with one or more genes in the ERN gene list are kept. Calls with no overlapping ERN genes are therefore removed. Two modes of filtering are possible: combined and singular ERN filtering. With combined ERN filtering, all ERN genelists are combined into one single list. This mode can be used by adding the `-c` flag to the command. With singular ERN filtering each sample is only the list appropriate to the sample is used to keep or remove calls.

__Required parameters__
* [-i / --infile]: Path to the combined CCRS calls file to filter
* [-e / --erndir]: Path to directory containing the ERN gene lists
* [-o / --outdir]: Directory path to write output files to
* [-op / --out-prefix]: Prefix to use for the output file(s)
* [-s / --samples-to-ern]: Table containing samples E-numbers linked to their ERN gene list

__Optional parameters__
* [-p / --percentage-overlap]: Required minimum percentage of overlap with an ERN gene (default=25). This parameter is currently not used.
* [-c / --combine-ern]: Whether to combine the multiple ERN gene lists into one list. (This is used for the Combined ERN filtering step)

__Usages__
```
python ern_filtering_2.py \
	-i /path/to/combined_ccrs/batch1.called.seg \
	-e /path/to/ern_genelists/ \
	-o /path/to/combined_ern_filtered_output/ \
	-op /batch1 \
	-s /path/to/samples_to_ern.txt
	-c
```

```
python ern_filtering.py \
	-i /path/to/combined_ccrs/batch1.called.seg \
	-e /path/to/ern_genelists/ \
	-o /path/to/singular_ern_filtered_output/ \
	-op /batch1 \
	-s /path/to/samples_to_ern.txt
```


## erngene.py
Contains the ErnGene class that `ern_filtering.py` uses to store ERN gene data.


## filter_ccrs_neutrals.py
Filters neutral 'calls' (segments with '0' that aren't deletions or duplications) and all calls from outlier samples from a set of CCRS calls file.

__Required parameters__
* [-i / --indir]: Directory containing CCRS call files (.igv.seg)
* [-o / --outdir]: Directory to write filtered CCRS call files to (.igv.seg)
* [-de / --deletion-outliers]: File containing deletion outlier sample names
* [-du / --duplication-outliers]: File containing duplication outlier sample names

__Usage__
```
python filter_ccrs_neutrals.py \
	-i /path/to/ccrs/batch1 \
	-o /path/to/filtered_neutrals/batch1/ \
	-de /path/to/batch1_del_outliers.txt \
	-du /path/to/batch1_dup_outliers.txt
```


## frequency_annotation.py
Determines the frequency of each call in two ways. First the occurrence and frequency of the call itself is determined. Second, the group occurrence and frequency is determined for each call. Groups of calls are created in several steps. First all calls are grouped by start position, and again by ending position. For groups of two or more calls the median call, determined by call length, is chosen as a representative. All calls overlapping at least x% or more with the representative call form an xpgroup. For these xpgroups, the occurrences and frequencies are determined by adding all individual occurrences. If a call is part of more than one xpgroup, it will receive the highest group occurrence and frequency. By default it will only annotate calls with the individual and group occurrence and frequency. By adding the `-f` flag, calls with a group occurrence and frequency that satifsy the set minimum number of calls and minimum percentage of samples will be removed.

Calls are annotated with the group label that gave the group occurrence and frequency.
There are six different group labels:
* ed: End position, decrease
* ei: End position, increase
* sd: Start position, decrease
* si: Start position, increase
* xd: Independent group, decrease
* xi: Independent group, increase

The labels are based on the group representative. If the xpgroup representative is a duplication call and originated from an end position group, the calls within the group will get the label 'ei' followed by a number (if that group has the highest frequency). If the representative is a deletion and originated from a start position group, the calls within the group will get the label 'sd' followed by a number.

__Required parameters__
* [-i / --infile]: Path to the combined CCRS file
* [-n / --numofsamples]: Number of processed samples for the BED file (batch)
* [-o / --outdir]: Path to the directory to write the new frequency annotated/filtered files to
* [-op / --outprefix]: Prefix to use for the output files

__Optional parameters__
* [-p / --percent-overlap]: Minimum percentage overlap required between calls and group representatives
* [-c / --minimum-calls]: Minimum number of calls required for calls from a group to be filtered (default=3).
* [-s / --minimum-samples]: Minimum percentage of samples required for calls from a group to be filtered (default=5.00).
* [-f / --filter-commonxp]: Remove 
* [-m / --maximum-size]: (default=5000000)


__Usage__
```
python frequency_annotation.py \
	-i /path/to/combined_ccrs/batch1.called \
	-n 178 \
	-o /path/to/frequency_output/ \
	-op batch1
```


## read_combined_ccrs.py
This script only contains the function to read a combined CCRS file, which is used in various other filtering and annotation scripts. It can read a combined CCRS file from different result processing steps, before or after columns such as CCRS_Frequency and Conrad_Frequency have been added. The method returns the header line and read data as CcrsCall objects saved per sample, per chromosome.


## solverdcall.py
Contains the SolveRdCall class used by compare_results.py to store basic data from other Solve-RD calls made by ClinCNV, Conifer and VarGenius.
