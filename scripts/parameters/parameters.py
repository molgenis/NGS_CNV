def get_classification_parameters(tool_choices):
    cmd_args = argparse.ArgumentParser()
    cmd_args.add_argument("-t", "--tool", type=str, dest="tool", choices=tool_choices, help="Name of the tool")
    cmd_args.add_argument("-s", "--samples", type=str, dest="samples", help="Path to samples file")
    cmd_args.add_argument("-e", "--exons", type=str, dest="exons", help="Path to exons file")
    cmd_args.add_argument("-p", "--probes", type=str, dest="probes", help="Path to probes file")
    cmd_args.add_argument("-i", "--indir", type=str, dest="indir", help="Path to input folder or file")
    cmd_args.add_argument("-a", "--array", type=str, dest="array", help="Path to array calls file")
    cmd_args.add_argument("-o", "--output", type=str, dest="output", help="Path to write output file to")
    cmd_args.add_argument("-op", "--outprefix", type=str, dest="outprefix", help="Prefix to use for output files")
    cmd_args.add_argument("-np", "--number-of-probes", type=int, dest="numofprobes", default=10, help="Required minimum number of probes to be Array Informative")
    cmd_args.add_argument("-ne", "--number-of-exons", type=int, dest="numofexons", default=3, help="Required minimum number of exons to be WES Informative")
    cmd_args.add_argument("-po", "--percent-overlap", type=int, dest="percentoverlap", default=50, help="Percentage overlap required between Array en WES CNVs")
    cmd_args.add_argument("-fn", "--filter-neutrals", dest="filterneutrals", action="store_true", help="Filter out neutrals?")
    cmd_args.add_argument("-cs", "--cnv-size", type=int, dest="cnvsize", help="Minimum GATK4 CNV size")
    return vars(cmd_args.parse_args())


def get_filter_parameters(tool_choices):
    filter_args = argparse.ArgumentParser()
    filter_args.add_argument("-t", "--tool", dest="tool", type=str, choices=tool_choices, help="Type of filtering to perform")
    filter_args.add_argument("-i", "--infile", dest="infile", type=str, help="Path to input file")
    filter_args.add_argument("-o", "--outfile", dest="outfile", type=str, help="Path to output file")
    filter_args.add_argument("-c", "--conradfile", dest="conradfile", type=str, help="Path to file with Conrad CNVs")
    filter_args.add_argumemt("-e", "--exonfile", dest="exonfile", type=str, help="Path to exon BED file")
    filter_args.add_argument("-u", "--commoncnvs", dest="commoncnvs", type=str, help="Path to file with Common CNVs")
    return vars(filter_args.parse_args())


def get_totals_parameters(tool_choices):
    totals_args = argparse.ArgumentParser()
    totals_args.add_argument("-t", "--tool", dest="tool", type=str, choices=tool_choices, required=True, help="Type of totals to collect")
    totals_args.add_argument("-i", "--infile", dest="infile", type=str, help="Path to input file")
    totals_args.add_argument("-a", "--arrayfile", dest="arrayfile", type=str, help="Path to file with array")
    totals_args.add_argument("-o", "--outfile", dest="outfile", type=str, help="")
    totals_args.add_argument("-op", "--outprefix", dest="outprefix", type=str, help="")
    return vars(totals_args.parse_args())


def get_util_parameters(tool_choices):
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


def required_parameters_set(tooltouse, filterparams, required_params):
    missing_parameters = []
    if tooltouse not in required_params:
        print("No valid tool selected")
        missing_parameters.append("tool")
    else:
        for paramname in required_params[tooltouse]:
            if filterparams[paramname] is None:
                missing_parameters.append(paramname)
    return missing_parameters


def parameters_are_ok(paramvalues, required_params, param_types):
    incorrect_parameters = []
    missing_parameters = required_params_set(paramvalues["tool"], paramvalues, required_params)

    if len(missing_parameters) == 0:
        for paramname in required_params[paramvalues["tool"]]:
            if param_types[paramname] == "inputfile":
                if not os.path.isfile(paramvalues[paramname]):
                    incorrect_parameters.append(paramname)
            elif param_types[paramname] == "outputfile":
                outputdir = "/".join(paramvalues[paramname].split("/")[0:-1])
                if not os.path.isdir(outputdir):
                    incorrect_parameters.append(paramname)
    else:
        incorrect_parameters.extend(missing_parameters)
    return incorrect_parameters
