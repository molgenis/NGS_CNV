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
    cmd_args.add_argument("-if", "infile", type=str, dest="infile", help="Path to input file")
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
    filter_args.add_argumemt("-e", "--exonfile", dest="exonfile", type=str, help="Path to exon BED file")
    filter_args.add_argument("-u", "--commoncnvs", dest="commoncnvs", type=str, help="Path to file with Common CNVs")
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
    totals_args.add_argument("-a", "--arrayfile", dest="arrayfile", type=str, help="Path to file with array")
    totals_args.add_argument("-o", "--outfile", dest="outfile", type=str, help="Path to write output to (file or dir)")
    totals_args.add_argument("-op", "--outprefix", dest="outprefix", type=str, help="Prefix to use for output file names")
    totals_args.add_argument("-po", "--percent-overlap", dest="percentoverlap", type=int, help="Minimum required percentage overlap")
    return vars(totals_args.parse_args())


def get_util_parameters(tool_choices):
    """Define and return

    Parameters
    ----------
    tool_choices : list of str
    """
    util_args = argparse.ArgumentParser()
    util_args.add_argument("-t", "--tool", dest="tool", type=str, choices=tool_choices, required=True, help="Specified util/tool to run")
    util_args.add_argument("-i", "--infile", dest="infile", type=str, help="Path to non interval or allelic input file")
    util_args.add_argument("-al", "--allelicfile", dest="allelicfile", type=str, help="Path to GATK4 allelic interval file")
    util_args.add_argument("-il", "--intervallist", dest="intervallist", type=str, help="Path to intervallist file")
    util_args.add_argument("-it", "--intype", dest="intype", type=str, required=False, choices=INTYPE_CHOICES, help="Input file type")
    util_args.add_argument("-o", "--outfile", dest="outfile", type=str, help="Path to output file")
    util_args.add_argument("-r", "--region", dest="region", nargs="+", type=str, help="Region(s) to use.")
    util_args.add_argument("-p", "--padding", dest="padding", default=10000, type=int, help="Amount of padding to add left and right of region(s)")
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
