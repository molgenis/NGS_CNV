import os
import sys
import argparse

# Import generate_totals scripts
import generate_totals.total_calls as gttc

# Import parameter script
import parameters.parameters as parpar

# Import util scripts
import utils.filereaders as ufr
import utils.filewriters as ufw


TOOL_CHOICES = ["arraycnv", "classification", "numofcalls", "numofna"]
REQUIRED_PARAMS = {"arraycnv": ["infile", "arrayfile"],
                   "classification": [],
                   "numofcalls": ["infile"],
                   "numofna": ["infile"]}
OPTIONAL_PARAMS = {}
PARAM_TYPES = {"infile": "inputfile",
               "arrayfile": "inputfile",
               "outfile": "outputfile",
               "outprefix": "string"}


def main():
    totalsparams = parpar.get_totals_params(TOOL_CHOICES)
    if parpar.parameters_are_ok(totalsparams, REQUIRED_PARAMS, PARAM_TYPES):
        if totalsparams["tool"] == "arraycnv":
            
        if totalsparams["tool"] == "classification":

        # Gather the total number of GATK4 Calls
        if totalsparams["tool"] == "numofcalls":
            gttc.get_total_calls(totalsparams["infile"])

        # Gather the total number of GATK4 calls that overlap with an array CNV
        if totalsparams["tool"] == "numofna":
            gttc.get_total_calls(totalsparams["infile"], True)


if __name__ == "__main__":
    main()
