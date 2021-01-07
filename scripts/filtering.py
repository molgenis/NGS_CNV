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
import filtering.naremoval as fnar
import filtering.size_filtering as fsf


# Create some general variables
TOOL_CHOICES = ["ccnvfiltering", "conradfiltering", "nafiltering", "sizefiltering"]
REQUIRED_PARAMS = {"ccnvfiltering": ["commoncnvs", "infile", "outfile"],
                   "conradfiltering": ["conradfile", "exonfile", "infile", "outfile"],
                   "nafiltering": ["infile", "outfile"],
                   "sizefiltering": ["infile", "outfile", "colname", "cnvsize"]}
OPTIONAL_PARAMS = {}
PARAM_TYPES = {"infile": "inputfile",
               "outfile": "outputfile",
               "conradfile",: "inputfile",}
               "exonfile",: "inputfile",}
               "commoncnvs",: "inputfile",
               "colname": "string",
               "cnvsize": "integer"}
TOOL_USAGE = {"ccnvfiltering": "python filtering.py -t ccnvfiltering -u commoncnvs.txt -i cnv_classifications.txt -o filtered_ccnv_classifications.txt",
              "conradfiltering": "python filtering.py -t conradfiltering -c cornad_cnvs.txt -e exons.bed -i cnv_classifications.txt -o filtered_cnv_classifications.txt",
              "nafitlering": "python filtering.py -t nafiltering -i cnv_classifications.txt -o nafiltered_cnv_classifications.txt",
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
            fnar.remove_nas(filterparams["infile"], filterparams["outfile"])

        # Filter CNV calls by size
        if filterparams["tool"] = "sizefiltering":
            run_sizefiltering(filterparams)
    else:
        print(f"Missing parameters: {incorrect_params}")
        parpar.display_tool_usage(filterparams["tool"], TOOL_USAGE)


def run_ccnvfiltering(filterparams):
    ccnvdata = ufr.read_umcg_common_cnv_file(filterparams["commoncnvs"])
    fuccf.filter_classification_results(filterparams["infile"], ccnvdata, filterparams["outfile"])


def run_conradfiltering(filterparams):
    gatkresultdata = fcf.read_classification_file(filterparams["infile"])
    conradcnvs = fcf.read_conrad_data(filterparams["conradfile"])
    fcf.add_qxte_exon_data(filterparams["exonfile"], conradcnvs)
    gatkconradcnv = fcf.determine_gatk_conrad_overlaps(gatkresultdata, conradcnvs)
    filteredgatk = fcf.determine_gatk_filtered_by_conrad(gatkconradcnv)
    filewritten = fcf.write_conrad_filtered_gatkcalls(gatkresults, filteredgatk, tg_params["outfile"])


def run_sizefiltering(filterparams):
    print("...Filtering CNV calls by size...")
    fsf.filter_by_size(filterparams["infile"], filterparams["outfile"], filterparams["colname"], filterparams["cnvsize"])


if __name__ == "__main__":
    main()
    print("DONE!")
