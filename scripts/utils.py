import os
import sys
import argparse

# Import util scripts
import utils.filereaders as ufr
import utils.filewriters as ufw
import utils.filter_xy_from_intervallist as ufxyfi
import utils.fix_array_cnvs as ufac
import utils.get_snp_log2_ratios as ugsl2r
import utils.select_plot_region as uspr


# Create some general variables
TOOL_CHOICES = ["filterxy", "fixarray", "getsnplog2ratios", "selectplotregion"]
REQUIRED_ARGS = {"filterxy": ["intervallist", "outfile"],
                 "fixarray": ["infile", "outfile"],
                 "getsnplog2ratios": ["allelicfile", "intervallist", "outfile"],
                 "selectplotregion": ["intervallist", "outfile", "padding", "region"]}
OPTIONAL_ARGS = {}
XY_CHROMS = ("X", "x", "Y", "y", "chrX", "chrx", "chrY", "chry")
INTYPE_CHOICES = ["cac", "tsv", "seg"]


def main():
    utilparams = get_params()
    if params_are_ok(utilparams):
        # Filter X and Y chromosomes from an interval list file
        if utilparams["tool"] == "filterxy":
            ufxyfi.filter_intervallist(utilparams["intervallist"], utilparams["outfile"])

        # Fix incorrect array cnv notation (current script versions do not need this!)
        if utilparams["tool"] == "fixarray":
            ufac.fix_array_cnv_notation(utilparams["infile"], utilparams["outfile"])

        # Selet data for a specified plot region
        if utilparams["tool"] == "selectplotregion":
            genomicregions = uspr.parse_regions(utilparams["region"], utilparams["padding"])
            intervaldata = uspr.extract_intervals(utilparams["infile"], genomicregions)
            
            # Check which type of output file to write (should correspond to the input file type)
            if utilparams["intype"] == "tsv":
                ufw.write_tsv_output_file(utilparams["outfile"], intervaldata)
            elif utilparams["intype"] == "seg":
                ufw.write_seg_output_file(utilparams["outfile"], intervaldata)
            elif utilparams["intype"] == "cac":
                ufw.write_cac_output_file(utilparams["outfile"], intervaldata)

        # Get interval log2 ratios for SNP positions
        if utilparams["tool"] == "getsnplog2ratios":
            allelicdata = ufr.read_allelic_data(utilparams["allelicfile"])
            intervaldata = ufr.read_interval_data(utilparams["intervallist"])
            ugsl2r.collect_allelic_log2_values(allelicdata, intervaldata, utilparams["outfile"])


def get_params():
    """Define and receive CLI set parameters.
    
    Returns
    -------
    dict
        Set parameter values
    """
    util_args = argparse.ArgumentParser()
    util_args.add_argument("-t", "--tool", dest="tool", type=str, choices=TOOL_CHOICES, required=True, help="Specified util/tool to run")
    util_args.add_argument("-i", "--infile", dest="infile", type=str, help="Path to non interval or allelic input file")
    util_args.add_argument("-al", "--allelicfile", dest="allelicfile", type=str, help="Path to GATK4 allelic interval file")
    util_args.add_argument("-il", "--intervallist", dest="intervallist", type=str, help="Path to intervallist file")
    util_args.add_argument("-it", "--intype", dest="intype", type=str, required=False, choices=INTYPE_CHOICES, help="Input file type")
    util_args.add_argument("-o", "--outfile", dest="outfile", type=str, help="Path to output file")
    util_args.add_argument("-r", "--region", dest="region", nargs="+", type=str, help="Region(s) to use.")
    util_args.add_argument("-p", "--padding", dest="padding", default=10000, type=int, help="Amount of padding to add left and right of region(s)")
    return vars(util_args.parse_args())


def required_params_set(tooltouse, utilparamvals):
    """Check and return whether the tool required parameters are set.

    Parameters
    ----------
    tooltouse : str
        Specific tool to run
    utilparamvals : dict
        Set parameter values

    Returns
    -------
    list of str
        Empty list if all parameters are set ; list with missing parameter names otherwise
    """
    missingparams = []
    for paramname in REQUIRED_ARGS[tooltouse]:
        if utilparamvals[paramval] is None:
            missingparams.append(paramname)
    return missingparams


def params_are_ok(argvalues):
    mising_parameters = required_params_set(argvalues["tool"], argvalues)
    if len(missing_parameters) > 0:
        print(f"Missing parameters {missing_parameters}")
        return False
    
    # Perform further parameter checking (such as file exists)
    
    return True


if __name__ == "__main__":
    main()
