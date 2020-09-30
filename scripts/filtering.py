import os
import sys
import argparse


# Import required classes
from classes.conradcnv import ConradCnv
from classes.conradexon import ConradExon
from classes.exon import Exon
from classes.gatkcall import GatkCall
from classes.umcgcommoncnv import UmcgCommonCnv

# Import util scripts
import utils.filereaders as ufr

# Import filtering scripts
import filtering.ccnv_filtering as fuccf
import filtering.conrad_filtering as fcf
import filtering.naremoval as fnar


# Create some general variables
TOOL_CHOICES = ["ccnvfiltering", "conradfiltering", "nafiltering"]
REQUIRED_PARAMS = {"ccnvfiltering": ["commoncnvs", "infile", "outfile"],
                   "conradfiltering": ["conradfile", "exonfile", "infile", "outfile"],
                   "nafiltering": ["infile", "outfile"]}
OPTIONAL_PARAMS = {}
PARAM_TYPES = {"infile": "inputfile",
               "outfile": "outputfile",
               "conradfile",: "inputfile",}
               "exonfile",: "inputfile",}
               "commoncnvs",: "inputfile"}


def main():
    filterparams = get_params()
    missingparams = required_params_set(filterparams["tool"], filterparams)
    if len(missingparams) == 0:
        # Perform filtering with a Common CNV list.
        if filterparams["tool"] == "ccnvfiltering":
            ccnvdata = ufr.read_umcg_common_cnv_file(filterparams["commoncnvs"])
            fuccf.filter_classification_results(filterparams["infile"], ccnvdata, filterparams["outfile"])

        # Perform filtering with Conrad CNVs
        if filterparams["tool"] == "conradfiltering":
            gatkresultdata = fcf.read_classification_file(filterparams["infile"])
            conradcnvs = fcf.read_conrad_data(filterparams["conradfile"])
            fcf.add_qxte_exon_data(filterparams["exonfile"], conradcnvs)
            gatkconradcnv = fcf.determine_gatk_conrad_overlaps(gatkresultdata, conradcnvs)
            filteredgatk = fcf.determine_gatk_filtered_by_conrad(gatkconradcnv)
            filewritten = fcf.write_conrad_filtered_gatkcalls(gatkresults, filteredgatk, tg_params["outfile"])

        # Remove GATK4 CNV calls without an overlapping Array CNV call
        if filterparams["tool"] == "nafiltering":
            fnar.remove_nas(filterparams["infile"], filterparams["outfile"])
    else:
        print(f"Missing parameters {missingparams}")


def get_params():
    """Define and return set CLI parameter values.

    Returns
    -------
    dict
        CLI set parameter values
    """
    filter_args = argparse.ArgumentParser()
    filter_args.add_argument("-t", "--tool", dest="tool", type=str, choices=TOOL_CHOICES, help="Type of filtering to perform")
    filter_args.add_argument("-i", "--infile", dest="infile", type=str, help="Path to input file")
    filter_args.add_argument("-o", "--outfile", dest="outfile", type=str, help="Path to output file")
    filter_args.add_argument("-c", "--conradfile", dest="conradfile", type=str, help="Path to file with Conrad CNVs")
    filter_args.add_argumemt("-e", "--exonfile", dest="exonfile", type=str, help="Path to exon BED file")
    filter_args.add_argument("-u", "--commoncnvs", dest="commoncnvs", type=str, help="Path to file with Common CNVs")
    return vars(filter_args.parse_args())


def required_params_set(tooltouse, filterparams):
    """Check if all required parameters are set.

    Parameters
    ----------
    tooltouse : str
        Specific filtering to perform
    filterparams : dict
        All set parameter values
    """
    missing_parameters = []
    if tooltouse not in REQUIRED_PARAMS:
        print("No valid tool selected")
        missing_parameters.append("tool")
    else:
        for paramname in REQUIRED_PARAMS[tooltouse]:
            if filterparams[paramname] is None:
                missing_parameters.append(paramname)
    return missing_parameters


def params_are_ok(paramvalues):
    """Check if parameters are set and are ok.

    Parameters
    ----------
    paramvalues : dict
        Set parameter values to check

    Returns
    -------
    incorrect_parameters : list of str
        Names of not set or incorrect parameters ; empty list if everything is ok
    """
    incorrect_parameters = []
    missing_parameters = required_params_set(paramvalues["tool"], paramvalues)

    if len(missing_parameters) == 0:
        for paramname in REQUIRED_PARAMS[paramvalues["tool"]]:
            if PARAM_TYPES[paramname] == "inputfile":
                if not os.path.isfile(paramvalues[paramname]):
                    incorrect_parameters.append(paramname)
            elif PARAM_TYPES[paramname] == "outputfile":
                outputdir = "/".join(paramvalues[paramname].split("/")[0:-1])
                if not os.path.isdir(outputdir):
                    incorrect_parameters.append(paramname)
    else:
        incorrect_parameters.extend(missing_parameters)
    return incorrect_parameters


if __name__ == "__main__":
    main()
