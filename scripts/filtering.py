import os
import sys
import argparse


# Import required classes
from classes.conradcnv import ConradCnv
from classes.conradexon import ConradExon
from classes.exon import Exon
from classes.gatkcall import GatkCall
from classes.umcgcommoncnv import UmcgCommonCnv

# Import parameters scripts
import parameters.parameters as parpar

# Import util scripts
import utils.filereaders as ufr

# Import filtering scripts
import filtering.ccnv_filtering as fuccf
import filtering.conrad_filtering as fcf
import filtering.na_filtering as fnaf
import filtering.naremoval as fnar
import filtering.size_filtering as fsf


# Create some general variables
TOOL_CHOICES = ["ccnvfiltering", "conradfiltering", "nafiltering", "sizefiltering"]
REQUIRED_PARAMS = {"ccnvfiltering": ["commoncnvs", "infile", "outfile"],
                   "conradfiltering": ["conradfile", "exonfile", "infile", "outfile"],
                   "nafiltering": ["infile", "outfile", "colname"],
                   "sizefiltering": ["infile", "outfile", "colname", "cnvsize"]}
OPTIONAL_PARAMS = {}
PARAM_TYPES = {"infile": "inputfile",
               "outfile": "outputfile",
               "conradfile": "inputfile",
               "exonfile": "inputfile",
               "commoncnvs": "inputfile",
               "colname": "string",
               "cnvsize": "integer"}
TOOL_USAGE = {"ccnvfiltering": "python filtering.py -t ccnvfiltering -u commoncnvs.txt -i cnv_classifications.txt -o filtered_ccnv_classifications.txt",
              "conradfiltering": "python filtering.py -t conradfiltering -c cornad_cnvs.txt -e exons.bed -i cnv_classifications.txt -o filtered_cnv_classifications.txt",
              "nafitlering": "python filtering.py -t nafiltering -i cnv_classifications.txt -o nafiltered_cnv_classifications.txt -C ArrayCNV",
              "sizefiltering": "python filtering.py -t sizefiltering -i cnv_classifications.txt -o filtered_classifications.txt -C gatk4Call -S 50000"}


def main():
    filterparams = parpar.get_filtering_parameters(TOOL_CHOICES)
    incorrect_params = parpar.parameters_are_ok(filterparams, REQUIRED_PARAMS, PARAM_TYPES)
    if len(incorrect_params) == 0:
        # Perform filtering with a Common CNV list.
        if filterparams["tool"] == "ccnvfiltering":
            run_ccnvfiltering(filterparams)

        # Perform filtering with Conrad CNVs
        if filterparams["tool"] == "conradfiltering":
            run_conradfiltering(filterparams)

        # Remove GATK4 CNV calls without an overlapping Array CNV call
        if filterparams["tool"] == "nafiltering":
            # fnar.remove_nas(filterparams["infile"], filterparams["outfile"])
            run_nafiltering(filterparams)

        # Filter CNV calls by size
        if filterparams["tool"] == "sizefiltering":
            run_sizefiltering(filterparams)
    else:
        print(f"Missing parameters: {incorrect_params}")
        parpar.display_tool_usage(filterparams["tool"], TOOL_USAGE)


def run_ccnvfiltering(filterparams):
    ccnvdata = ufr.read_umcg_common_cnv_file(filterparams["commoncnvs"])
    fuccf.filter_classification_results(filterparams["infile"], ccnvdata, filterparams["outfile"])


def run_conradfiltering(filterparams):
    print("...Reading classification data...")
    gatkresultdata = fcf.read_classification_file(filterparams["infile"])
    print("...Reading Conrad CNV data...")
    conradcnvs = fcf.read_conrad_data(filterparams["conradfile"])
    print("...Linking exons fropm BED to Conrad CNVs...")
    fcf.add_qxte_exon_data(filterparams["exonfile"], conradcnvs)
    print("...Determining classification calls Conrad overlap...")
    gatkconradcnv = fcf.determine_gatk_conrad_overlaps(gatkresultdata, conradcnvs)
    print("...Filtering classification calls with Conrad CNVs...")
    filteredgatk = fcf.determine_gatk_filtered_by_conrad(gatkconradcnv)
    filewritten = fcf.write_conrad_filtered_gatkcalls(gatkresultdata, filteredgatk, filterparams["outfile"])
    print(f"...Wrote output file: {filewritten}...")


def run_nafiltering(filterparams):
    print("...Filtering NAs...")
    file_written = fnaf.filter_na(filterparams["infile"], filterparams["outfile"], filterparams["colname"])
    print(f"...Wrote output file: {file_written}...")


def run_sizefiltering(filterparams):
    print("...Filtering CNV calls by size...")
    file_written = fsf.filter_by_size(filterparams["infile"], filterparams["outfile"], filterparams["colname"], filterparams["cnvsize"])
    print(f"...Wrote output file: {file_written}...")


if __name__ == "__main__":
    main()
    print("DONE!")
