import os
import sys
import argparse

# Import parameter scripts
import parameters.parameters as parpar

# Import util scripts
import utils.arraycnv_bedregions as acbr
import utils.filereaders as ufr
import utils.filewriters as ufw
import utils.filter_xy_from_intervallist as ufxyfi
import utils.fix_array_cnvs as ufac
import utils.get_snp_log2_ratios as ugsl2r
import utils.select_plot_region as uspr


# Create some general variables
TOOL_CHOICES = ["arraybedregion", "filterxy", "fixarray", "getsnplog2ratios", "selectplotregion"]
REQUIRED_PARAMS = {"arraybedregion": ["bedfile", "region"],
                   "filterxy": ["intervallist", "outfile"],
                   "fixarray": ["infile", "outfile"],
                   "getsnplog2ratios": ["allelicfile", "intervallist", "outfile"],
                   "selectplotregion": ["intervallist", "outfile", "padding", "region"]}
OPTIONAL_PARAMS = {}
PARAM_TYPES = {"infile": "inputfile",
               "allelicfile": "inputfile",
               "intervallist": "inputfile",
               "intype": "string",
               "outfile": "outputfile",
               "region": "string",
               "padding": "integer",
               "bedfile": "inputfile"}
XY_CHROMS = ("X", "x", "Y", "y", "chrX", "chrx", "chrY", "chry")
INTYPE_CHOICES = ["cac", "tsv", "seg"]
TOOL_USAGE = {"filterxy": "python utils.py -il intervallist.txt -o no_xy_intervallist.txt",
              "fixarray": "python utils.py -i cnv_classifications.txt -o fixedcnv_classifications.txt",
              "getsnplog2ratios": "python utils.py -al allelicratios.csv -il intervallist.txt -o snplog2ratios.txt",
              "selectplotregion": "python utils.py -il intervallist.txt -o region_interval_data.txt -p 1000 -r chr1:100-1000",}


def main():
    utilparams = parpar.get_util_parameters(TOOL_CHOICES)
    incorrectparams = parpar.parameters_are_ok(utilparams, REQUIRED_PARAMS, PARAM_TYPES)

    if len(incorrectparams) == 0:
        # Display BED regions overlapping with an array CNV region
        if utilparams["tool"] == "arraybedregion":
            run_arraycnv_bedregions(utilparams)

        # Filter X and Y chromosomes from an interval list file
        if utilparams["tool"] == "filterxy":
            ufxyfi.filter_intervallist(utilparams["intervallist"], utilparams["outfile"])

        # Fix incorrect array cnv notation (current script versions do not need this!)
        if utilparams["tool"] == "fixarray":
            ufac.fix_array_cnv_notation(utilparams["infile"], utilparams["outfile"])

        # Selet data for a specified plot region
        if utilparams["tool"] == "selectplotregion":
            run_selectplotregion(utilparams)

        # Get interval log2 ratios for SNP positions
        if utilparams["tool"] == "getsnplog2ratios":
            run_getsnplog2ratios(utilparams)
    else:
        print(f"The following parameters are incorrect: {incorrectparams}")
        parpar.display_tool_usage(utilparams["tool"], TOOL_USAGE)


def run_selectplotregion(utilparams):
    genomicregions = uspr.parse_regions(utilparams["region"], utilparams["padding"])
    intervaldata = uspr.extract_intervals(utilparams["infile"], genomicregions)
    
    # Check which type of output file to write (should correspond to the input file type)
    if utilparams["intype"] == "tsv":
        ufw.write_tsv_output_file(utilparams["outfile"], intervaldata)
    elif utilparams["intype"] == "seg":
        ufw.write_seg_output_file(utilparams["outfile"], intervaldata)
    elif utilparams["intype"] == "cac":
        ufw.write_cac_output_file(utilparams["outfile"], intervaldata)


def run_getsnplog2ratios(utilparams):
    allelicdata = ufr.read_allelic_data(utilparams["allelicfile"])
    intervaldata = ufr.read_interval_data(utilparams["intervallist"])
    ugsl2r.collect_allelic_log2_values(allelicdata, intervaldata, utilparams["outfile"])


def run_arraycnv_bedregions(utilparams):
    print("...Reading BED file data...")
    bedfiledata = ufr.read_exon_data(utilparams["bedfile"])
    print("...Collecting overlapping BED regions...")
    arraybedregions = acbr.get_arraycnv_bedregions(bedfiledata, utilparams["region"])
    acbr.display_bedregions(arraybedregions, utilparams["region"])


if __name__ == "__main__":
    main()
