import argparse
import os

def get_classification_parameters(tool_choices):
    """Define and return set CLI parameters for the classification script.

    Parameters
    ----------
    tool_choices : list of str
        List of tools that can be chosen

    Returns
    -------
    dict
        Parameter values
    """
    cmd_args = argparse.ArgumentParser()
    cmd_args.add_argument("-t", "--tool", type=str, dest="tool", choices=tool_choices, help="Name of the tool to use")
    cmd_args.add_argument("-s", "--samples", type=str, dest="samples", help="Path to samples file")
    cmd_args.add_argument("-e", "--exonsfile", type=str, dest="exonsfile", help="Path to exons file")
    cmd_args.add_argument("-p", "--probesfile", type=str, dest="probesfile", help="Path to probes file")
    cmd_args.add_argument("-i", "--indir", type=str, dest="indir", help="Path to input folder")
    cmd_args.add_argument("-if", "--infile", type=str, dest="infile", help="Path to input file")
    cmd_args.add_argument("-a", "--arrayfile", type=str, dest="arrayfile", help="Path to array CNV calls file")
    cmd_args.add_argument("-o", "--output", type=str, dest="output", help="Path to folder to write output file to")
    cmd_args.add_argument("-op", "--outprefix", type=str, dest="outprefix", help="Prefix to use for output files")
    cmd_args.add_argument("-np", "--number-of-probes", type=int, dest="numofprobes", default=10, help="Minimum required number of probes to be Array Informative")
    cmd_args.add_argument("-ne", "--number-of-exons", type=int, dest="numofexons", default=3, help="Minimum required number of exons to be WES Informative")
    cmd_args.add_argument("-po", "--percent-overlap", type=int, dest="percentoverlap", default=50, help="Percentage overlap required between Array en WES CNVs to be considered True Positive")
    cmd_args.add_argument("-fn", "--filter-neutrals", dest="filterneutrals", action="store_true", help="Filter out neutrals calls (GATK4 specific)?")
    cmd_args.add_argument("-cs", "--cnv-size", type=int, dest="cnvsize", help="Minimum required CNV size")
    cmd_args.add_argument("-ed", "--edsamples", type=str, dest="exomedepthsamples", help="Path to table linking samples and ExomeDepth output files")
    return vars(cmd_args.parse_args())


def get_comparison_parameters(tool_choices):
    """Define and return set CLI parameters for the comparison script.

    Parameters
    ----------
    tool_choices : list of str
        List of comparisons that can be chosen

    Returns
    -------
    dict
        Parameter values for comparison script
    """
    compare_args = argparse.ArgumentParser()
    compare_args.add_argument("-t", "--tool", type=str, required=True, choices=tool_choices, dest="tool", help="Type of comparison to perform")
    compare_args.add_argument("-a", "--arraycnvs", type=str, dest="arraycnvs", help="Path to file with array CNVs")
    compare_args.add_argument("-c", "--conifer-file", type=str, dest="conifer-file", help="Path to Conifer classification file")
    compare_args.add_argument("-c2", "--conifer-file2", type=str, dest="conifer-file2", help="Path to second Conifer classification file")
    compare_args.add_argument("-e", "--exomedepth-file", type=str, dest="exomedepth-file", help="Path to ExomeDepth classification file")
    compare_args.add_argument("-e2", "--exomedepth-file2", type=str, dest="exomedepth-file2", help="Path to second ExomeDepth classification file")
    compare_args.add_argument("-g", "--gatk4-file", type=str, dest="gatk4-file", help="Path to GATK4 classification file")
    compare_args.add_argument("-g2", "--gatk4-file2", type=str, dest="gatk4-file2", help="Path to the second GATK4 classification file")
    compare_args.add_argument("-o", "--outdir", type=str, dest="outdir", help="Path to write the comparison output files to")
    compare_args.add_argument("-op", "--output-prefix", type=str, dest="output-prefix", help="Prefix to use for the output files")
    compare_args.add_argument("-s", "--sample-file", type=str, dest="sample-file", help="Path to sample table")
    compare_args.add_argument("-1", "--file1", type=str, dest="file1", help="Path to first CNV calling classification file")
    compare_args.add_argument("-2", "--file2", type=str, dest="file2", help="Path to second CNV calling classification file")
    compare_args.add_argument("-l1", "--label1", type=str, dest="label1", help="Label to use for the first tool")
    compare_args.add_argument("-l2", "--label2", type=str, dest="label2", help="Label to use for the second tool")
    compare_args.add_argument("--tp-per-acnv", dest="tp-per-acnv", action="store_true", help="Count TPs only per array CNV")
    return vars(compare_args.parse_args())


def get_dualbed_parameters(tool_choices):
    """Define and return set CLI parameters for the dualbed script.

    Parameters
    ----------
    tool_choices : list of str
        List of actions for dualBED that can be chosen
    """
    dualbed_args = argparse.ArgumentParser()
    dualbed_args.add_argument("-t", "--tool", type=str, dest="tool", choices=tool_choices, help="Type of action to perform")
    dualbed_args.add_argument("-1", "--infile1", type=str, dest="infile1", help="Path to first input file")
    dualbed_args.add_argument("-2", "--infile2", type=str, dest="infile2", help="Path to second input file")
    dualbed_args.add_argument("-o", "--outfile", type=str, dest="outfile", help="Path to write output file to")
    dualbed_args.add_argument("-od", "--outdir", type=str, dest="outdir", help="Path to output directory to write output files to")
    dualbed_args.add_argument("-op", "--output-prefix", type=str, dest="output-prefix", help="Output prefix to use")
    dualbed_args.add_argument("-po", "--percent-overlap", type=int, dest="percent-overlap", default=75, help="Minimal required percentage overlap")
    dualbed_args.add_argument("-nu", "--no-unique", dest="no-unique", action="store_true", help="Exclude calls labelled Unique")
    return vars(dualbed_args.parse_args())


def get_filtering_parameters(tool_choices):
    """Define and return set CLI parameters for the filtering script.

    Parameters
    ----------
    tool_choices : list of str
        List of tools that can be chosen

    Returns
    -------
    dict
        Parameter values for filter script
    """
    filter_args = argparse.ArgumentParser()
    filter_args.add_argument("-t", "--tool", dest="tool", type=str, choices=tool_choices, help="Type of filtering to perform")
    filter_args.add_argument("-i", "--infile", dest="infile", type=str, help="Path to input file")
    filter_args.add_argument("-o", "--outfile", dest="outfile", type=str, help="Path to output file")
    filter_args.add_argument("-c", "--conradfile", dest="conradfile", type=str, help="Path to file with Conrad CNVs")
    filter_args.add_argument("-e", "--exonfile", dest="exonfile", type=str, help="Path to exon BED file")
    filter_args.add_argument("-u", "--commoncnvs", dest="commoncnvs", type=str, help="Path to file with Common CNVs")
    filter_args.add_argument("-C", "--colname", dest="colname", type=str, help="Column name for size filtering")
    filter_args.add_argument("-S", "--cnvsize", dest="cnvsize", type=int, help="Minimum size for CNVs to be retained")
    return vars(filter_args.parse_args())


def get_totals_parameters(tool_choices):
    """Define and return set CLI parameters for the totals script.

    Parameters
    ----------
    tool_choices : list of str
        List of tools that can be chosen

    Returns
    -------
    dict
        Set parameter values for the totals script
    """
    totals_args = argparse.ArgumentParser()
    totals_args.add_argument("-t", "--tool", dest="tool", type=str, choices=tool_choices, required=True, help="Type of totals to collect")
    totals_args.add_argument("-i", "--infile", dest="infile", type=str, help="Path to input file")
    totals_args.add_argument("-2", "--infile2", dest="infile2", type=str, help="Patht op second input file")
    totals_args.add_argument("-a", "--arrayfile", dest="arrayfile", type=str, help="Path to file with array")
    totals_args.add_argument("-o", "--outfile", dest="outfile", type=str, help="Path to write output to (file or dir)")
    totals_args.add_argument("-op", "--outprefix", dest="outprefix", type=str, help="Prefix to use for output file names")
    totals_args.add_argument("-po", "--percent-overlap", dest="percentoverlap", type=int, help="Minimum required percentage overlap")
    totals_args.add_argument("--tp-per-acnv", dest="tp-per-acnv", action="store_true", help="Count TPs only per array CNV")
    return vars(totals_args.parse_args())


def get_util_parameters(tool_choices):
    """Define and return

    Parameters
    ----------
    tool_choices : list of str
    """
    INTYPE_CHOICES = ["cac", "tsv", "seg"]
    util_args = argparse.ArgumentParser()
    util_args.add_argument("-t", "--tool", dest="tool", type=str, choices=tool_choices, required=True, help="Specified util/tool to run")
    util_args.add_argument("-i", "--infile", dest="infile", type=str, help="Path to non interval or allelic input file")
    util_args.add_argument("-al", "--allelicfile", dest="allelicfile", type=str, help="Path to GATK4 allelic interval file")
    util_args.add_argument("-il", "--intervallist", dest="intervallist", type=str, help="Path to intervallist file")
    util_args.add_argument("-it", "--intype", dest="intype", type=str, required=False, choices=INTYPE_CHOICES, help="Input file type")
    util_args.add_argument("-o", "--outfile", dest="outfile", type=str, help="Path to output file")
    util_args.add_argument("-r", "--region", dest="region", type=str, help="Region(s) to use.")
    util_args.add_argument("-p", "--padding", dest="padding", default=10000, type=int, help="Amount of padding to add left and right of region(s)")
    util_args.add_argument("-b", "--bedfile", dest="bedfile", type=str, help="Path to BED file")
    return vars(util_args.parse_args())


def required_parameters_set(tooltouse, paramvalues, required_params):
    """Check and return whether tool required parameters have been set.

    Parameters
    ----------
    tooltouse : str
        Tool to check parameters for
    paramvalues : dict
        Set parameter values
    required_params : dict
        Required parameters for tools

    Returns
    -------
    list of str
        Names of missing parameters ; empty list if all required parameters have been set
    """
    missing_parameters = []
    if tooltouse not in required_params:
        print("No valid tool selected")
        missing_parameters.append("tool")
    else:
        for paramname in required_params[tooltouse]:
            if paramvalues[paramname] is None:
                missing_parameters.append(paramname)
    return missing_parameters


def parameters_are_ok(paramvalues, required_params, param_types):
    """Check parameter values by determining whether input files exist, etc.

    Parameters
    ----------
    paramvalues : dict
       Set parameter values 
    required_params : dict
    param_types : dict

    Returns
    -------
    incorrect_parameters : list of str
        Names of incorrect parameters ; empty list if all required parameters are ok.
    """
    incorrect_parameters = []
    missing_parameters = required_parameters_set(paramvalues["tool"], paramvalues, required_params)

    if len(missing_parameters) == 0:
        for paramname in required_params[paramvalues["tool"]]:
            
            # Check if a provided file exists.
            if param_types[paramname] == "inputfile":
                if not os.path.isfile(paramvalues[paramname]):
                    incorrect_parameters.append(paramname)

            # Check whether the folder to write an output file to exists.
            elif param_types[paramname] == "outputfile":
                outputdir = ""
                if '/' in paramvalues[paramname]:
                    outputdir = "/".join(paramvalues[paramname].split("/")[0:-1])
                elif '\\' in paramvalues[paramname]:
                    outputdir = "\\".join(paramvalues[paramname].split("\\")[0:-1])
                if not os.path.isdir(outputdir):
                    incorrect_parameters.append(paramname)

            # Check to make sure that the parameter that should be a string is a indeed a string and also not empty
            elif param_types[paramname] == "string":
                if not type(paramvalues[paramname]) == str:
                    incorrect_parameters.append(paramname)
                elif paramvalues[paramname] == "":
                    incorrect_parameters.append(paramname)

            # Check to make sure that the parameter that should be an integer is indeed an integer
            elif param_types[paramname] == "integer":
                if not type(paramvalues[paramname]) == int:
                    incorrect_parameters.append(paramname)
    else:
        incorrect_parameters.extend(missing_parameters)
    return incorrect_parameters


def display_tool_usage(tooltouse, tools_usage):
    """Display the usage command for a specific tool.

    Parameters
    ----------
    tooltouse : str
        Specific tool to display usage command for
    tools_usage : dict
        Example usage command for each tool
    """
    if tooltouse in tools_usage:
        print(f"{tooltouse} usage: {tools_usage[tooltouse]}")
    else:
        print("Selected tool does not exist :(")
