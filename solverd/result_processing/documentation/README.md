# Result processing documentation scripts

## calls_per_ern.py
Counts the number of CNV calls per sample for a single ERN. The output file craeted by `combine_calls_per_ern.py` us expected as input. The output file contains two columns: samplename and call count.

__Required parameter__
* [-i / --infile]: Path to the file (created with `combine_calls_per_ern.py`) containing all GATK4 CNV calls for a single ERN
* [-s / --samples-to-ern]: Path to the RD3 file linking sample names to ERNs
* [-o / --outfile]: Path ot write the output tabel to
* [-l / --label]: Specific ERN to count calls per sample for (choices: ["ERN EURO-NMD", "ERN GENTURIS", "ERN ITHACA", "ERN-RND", "UDN-Spain"])


## ccrs_num_of_calls.py


## ccrs_probes_per_call.py
Collects the number of probes from of each GATK4 CNV call. A combined CCRS file, created by `combine_ccrs.py` is expected as input. The number of probes are collected and written to files for duplications and deletions separately.

__Required parameters__
* [-i / --infile]: Path to combined CCRS file for a batch (BED)
* [-o / --outdir]: Directory to write output files to
* [-p / --prefix]: Prefix to use for output files
