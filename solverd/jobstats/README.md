# JobStats
Contains two scripts that can be used to get the average, median and stdev runtime of multiple slurm jobs, such as multiple CollectReadCounts or CollectAllelicCount jobs.

## 1: obtain_job_data.py
Is used to combine job information (JobID, Elapsed, AllocCPUS, AveCPU, ReqMem, MaxVMSize, MaxRSS, MaxDiskRead, MaxDiskWrite) from multiple slurm job .out files into a single tab separated file.

__Required parameters__
* [-d / --indir]: Path to directory containing multiple slurm jobs .out files
* [-o / --outfile]: Path to write combined job data file to.

__Usage__
```
python obtain_job_data.py \
	-d /path/to/crc_jobout_dir/ \
	-o /path/to/combined_crc_jobouts.txt
```


## 2: jobdatastats.py
Is used to calculate the average, median and stdev runtime for a job in seconds, minutes and hours.

__Required parameters__
* [-i / --infile]: Path to the combined jobout file created with `obtain_job_data.py`
* [-j / --jobname]: Label to use for when displaying the avg, med and stdev runtimes

__Usage__
```
python jobdatastats.py\
	-i /path/to/combined_crc_jobouts.txt
	-j batch2_crc
```
