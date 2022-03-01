#!/usr/bin/env python
import sys
from collections import Counter

# Import classes
from classes.arraycnv import ArrayCnv
from classes.fpoverlap import FpOverlap
from classes.fpregion import FpRegion
from classes.fpresult import FpResult
from classes.fpsummary import FpSummary

# Import generate_totals scripts
import generate_totals.array_cnvs as gtac
import generate_totals.classifications as gtcf
import generate_totals.fp_regions as gtfr
import generate_totals.total_calls as gttc
import generate_totals.dualbed_acnvs as gtdbac
import generate_totals.dualbed_ratios as gtdbr

# Import parameter script
import parameters.parameters as parpar

# Import shared methods scripts
import shared_methods.shared_methods as smsm

# Import util scripts
import utils.filereaders as ufr
import utils.filewriters as ufw


TOOL_CHOICES = ["arraycnv", "classification", "fpregions", "numofcalls", "numofna", "dualbed_arraycnv", "dualbed_classification"]
REQUIRED_PARAMS = {"arraycnv": ["arrayfile", "infile", "outfile", "outprefix"],
                   "classification": ["infile", "outfile"],
                   "fpregions": ["infile", "outfile", "outprefix"],
                   "numofcalls": ["infile"],
                   "numofna": ["infile"],
                   "dualbed_arraycnv": ["arrayfile", "infile", "outfile", "outprefix"],
                   "dualbed_classification": ["infile", "infile2", "outfile"]}
OPTIONAL_PARAMS = {}
PARAM_TYPES = {"infile": "inputfile",
               "infile2": "inputfile",
               "arrayfile": "inputfile",
               "outfile": "outputfile",
               "outprefix": "string",
               "percentoverlap": "integer"}
TOOL_USAGE = {"arraycnv": "python totals.py -t arraycnv -i cnv_classifications.txt -a arraycnvdata.txt",
              "classification": "python totals.py -t classification -i cnv_classifications.txt -o /path/to/outdir -op prefixname",
              "numofcalls": "python totals.py -t numofcalls -i cnv_classifications.txt",
              "numofna": "python totals.py -t numofna -i cnv_classifications.txt"}


def main():
    totalsparams = parpar.get_totals_parameters(TOOL_CHOICES)
    incorrect_parameters = parpar.parameters_are_ok(totalsparams, REQUIRED_PARAMS, PARAM_TYPES)
    if len(incorrect_parameters) == 0:
        # Determine the total number of array CNVs found and missed
        if totalsparams["tool"] == "arraycnv":
            run_arraycnv(totalsparams)

        # Gather the totals per classification label
        if totalsparams["tool"] == "classification":
            run_classification(totalsparams)

        # Gather the number of unique and duplicate FP regions
        if totalsparams["tool"] == "fpregions":
            run_fpregions(totalsparams)

        # Gather the total number of GATK4 Calls
        if totalsparams["tool"] == "numofcalls":
            gttc.get_total_calls(totalsparams["infile"])

        # Gather the total number of GATK4 calls that overlap with an array CNV
        if totalsparams["tool"] == "numofna":
            gttc.get_total_calls(totalsparams["infile"], True)

        # Gather the total number of array CNVs found and missed for dualBED
        if totalsparams["tool"] == "dualbed_arraycnv":
            run_dualbed_arraycnv(totalsparams)

        if totalsparams["tool"] == "dualbed_classification":
            run_dualbed_classification_ratios(totalsparams)
    else:
        print(f"Missing the following parameters: {incorrect_parameters}")


def run_arraycnv(totalsparams):
    """Determine the number of found and missed array CNVs.

    Parameters
    ----------
    totalsparams : dict
        Set CLI parameter values
    """
    outpath = totalsparams["outfile"] + "/" + totalsparams["outprefix"]
    arraydata = ufr.read_array_cnvs(totalsparams["arrayfile"])

    found_arraycnvs = gtac.determine_arraycnvs_found(totalsparams["infile"], arraydata)
    found_summary = gtac.summarize_arraycnv_types(found_arraycnvs)
    missed_arraycnvs = gtac.determine_arraycnvs_missed(found_arraycnvs, arraydata)
    missed_summary = gtac.summarize_arraycnv_types(missed_arraycnvs)

    ufw.write_missedfound_arraycnvs(found_arraycnvs, f"{outpath}_found_arraycnvs.txt")
    ufw.write_missedfound_summary(found_summary, f"{outpath}_found_arraycnvs_summary.txt")
    ufw.write_missedfound_arraycnvs(missed_arraycnvs, f"{outpath}_missed_arraycnvs.txt")
    ufw.write_missedfound_summary(missed_summary, f"{outpath}_missed_arraycnvs_summary.txt")


def run_dualbed_arraycnv(totalsparams):
    """Determine the number of found and missed array CNVs for dualBED results.

    Parameters
    ----------
    totalsparams : dict
        Set CLI parameters values
    """
    outpath = totalsparams["outfile"] + "/" + totalsparams["outprefix"]
    arraydata = ufr.read_array_cnvs(totalsparams["arrayfile"])

    # Determine the Shared, Overlapping and Unique array CNVs found
    s_found_acnvs = gtdbac.determine_dualbed_acnvs_found(totalsparams["infile"], arraydata, "Shared")
    shared_filter = gtdbac.make_sou_found_filter(s_found_acnvs)
    o_found_acnvs = gtdbac.determine_dualbed_acnvs_found_2(totalsparams["infile"], arraydata, "Overlapping", shared_filter)
    u_found_acnvs = gtdbac.determine_dualbed_acnvs_found(totalsparams["infile"], arraydata, "Unique")
    sou_found_summary = gtdbac.summarize_found_arraycnv_types(s_found_acnvs, o_found_acnvs, u_found_acnvs)

    # Determine the missed array CNVs
    missed_acnvs = gtdbac.determine_dualbed_acnvs_missed(arraydata, s_found_acnvs, o_found_acnvs, u_found_acnvs)
    missed_summary = gtdbac.summarize_missed_arraycnv_types(missed_acnvs)

    # Write the found and missed array CNV data to output files
    ufw.write_missedfound_arraycnvs(s_found_acnvs, f"{outpath}_shared_found_acnvs.txt")
    ufw.write_missedfound_arraycnvs(o_found_acnvs, f"{outpath}_overlapping_found_acnvs.txt")
    ufw.write_missedfound_arraycnvs(u_found_acnvs, f"{outpath}_unique_found_acnvs.txt")
    ufw.write_missedfound_arraycnvs(missed_acnvs, f"{outpath}_missed_acnvs.txt")

    # Write the found and missed summaries to output files
    gtdbac.write_found_summary(sou_found_summary, f"{outpath}_found_acnvs_summary.txt")
    gtdbac.write_missed_summary(missed_summary, f"{outpath}_missed_acnvs_summary.txt")


def run_classification(totalsparams):
    """Gather totals for classification labels.

    Parameters
    ----------
    totalsparams : dict
        Set CLI parameter values
    """
    print("...Reading classified CNV calls...")
    gatkresults = ufr.read_classification_file(totalsparams["infile"])
    print("...Generating classification label totals...")
    totalsdata = gtcf.generate_classification_totals(gatkresults, totalsparams["tp-per-acnv"])
    filewritten = ufw.write_classification_totals(totalsdata, totalsparams["outfile"])
    print(f"...Wrote outfile?: {filewritten}...")


def run_dualbed_classification_ratios(totalsparams):
    """Gather totals for dualBED classification ratios.

    Parameters
    ----------
    totalsparams : dict
        Set CLI parameter values
    """
    normal_data = gtdbr.read_dualbed_data(totalsparams["infile"])
    hc_data = gtdbr.read_dualbed_data(totalsparams["infile2"])

    # Generate totals for Shared dualBED calls
    shared_totals = gtcf.generate_classification_totals(normal_data["Shared"], totalsparams["tp-per-acnv"])

    # Generate totals for Overlapping dualBED calls
    overlapping_calls = gtdbr.merge_overlapping(normal_data["Overlapping"], hc_data["Overlapping"])
    sharedfilter = gtdbr.form_shared_filter(normal_data["Shared"])
    overlapping_totals = gtdbr.generate_overlapping_totals(overlapping_calls, sharedfilter, totalsparams["tp-per-acnv"])

    # Generate totals for Unique dualBED calls
    nunique_totals = gtcf.generate_classification_totals(normal_data["Unique"], totalsparams["tp-per-acnv"])
    hcunique_totals = gtcf.generate_classification_totals(hc_data["Unique"], totalsparams["tp-per-acnv"])

    # Write the results to file
    gtdbr.write_dualbed_ratios(totalsparams["outfile"], shared_totals, overlapping_totals, nunique_totals, hcunique_totals)


def run_fpregions(totalsparams):
    """Gather the total for False Positive regions.

    Parameters
    ----------
    totalsparams : dict
        Set CLI parameter values
    """
    # Read False Positive regions and count there occurences
    fpregions = ufr.read_fp_classifications(totalsparams["infile"])
    fpregion_counts = dict(Counter(fpregions))

    # Split the fpregion_counts into unique and duplicated regions. Also similar regions if needed
    dupfpregions = gtfr.get_duplicated_regions(fpregion_counts)
    unifpregions = gtfr.get_unique_regions(fpregion_counts)
    simfpregions = {}
    if totalsparams["percentoverlap"] is not None:
        if 0 < totalsparams["percentoverlap"] < 100:
            simfpregions = gtfr.determine_similar_regions(dupfpregions, unifpregions, totalsparams["percentoverlap"])

    # Write the output files with duplicate, similar and unique FP regions
    outpath = totalsparams["outfile"] + "/" + totalsparams["outprefix"]
    ufw.write_duplicate_fpregions(f"{outpath}_duplicate_fps.txt", dupfpregions)
    if len(simfpregions) > 0:
        ufw.write_similar_fpregions(f"{outpath}_similar_fps.txt", simfpregions)
        simfilterlist = gtfr.get_similar_filterlist(simfpregions)
        print(f"Found {len(simfilterlist)} similar False Positive regions")
        unifpregions = gtfr.filter_uniquefps_with_similars(simfilterlist, unifpregions)
    ufw.write_unique_fpregions(f"{outpath}_unique_fps.txt", unifpregions)


if __name__ == "__main__":
    main()
