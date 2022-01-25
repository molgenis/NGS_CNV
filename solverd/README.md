# Script folders:
* cnv_calling
* hcbedfile
* jobstats
* preprocessing
* result_processing


## cnv_calling
Contains script to generate GATK4 CNV calling jobs for each step (CollectReadCounts, CollectAllelicCounts, etc)


## hcbedfile
Contains a few scripts to hopfully help create the HC BDF file. These scripts have not been tested and may therefore not work yet. See more info in the README.


## jobstats
Contains two scripts to collect slurmjob runtimes from a folder containing the .out files. The obtain_job_data.py script can be used to collect all slurmjob information into a single file. The `jobdatastats.py` script can then be used to get the average, median and stdev job runtimes in seconds, minutes and hours.


## preprocessing
Contains scripts used to organize the data prior to CNV calling.


## result_processing
Contains multiple scripts to process the resulting CNV calls from the CCRS steps. Result processing include filtering away neutral segments, combining individual CCRS files into one and filtering CNV calls by overlap with ERN genes.
