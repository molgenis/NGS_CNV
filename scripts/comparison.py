import os
import argparse

# Import required classes
from classes.arraycnv import ArrayCnv
from classes.conifercall import ConiferCall
from classes.gatkcall import GatkCall

# Import comparison scripts
import comparison.comparison as comcom

# Import parameter scripts
import parameters.parameters as parpar

# Import util scripts
import utils.filereaders as ufr
import utils.filewriters as ufw


#Make some parameter defining variables
TOOL_CHOICES = ["gatk4_conifer", "gatk4_exomedepth", "conifer_exomedepth"]
REQUIRED_PARAMS = {"gatk4_conifer": ["gatk4-file", "conifer-file", "outdir", "output-prefix"],
                   "gatk4_exomedepth": ["gatk4-file", "exomedepth-file", "outdir", "output-prefix"],
                   "conifer_exomedepth": ["conifer-file", "exomedepth-file", "outdir", "output-prefix"]}
OPTIONAL_PARAMS = {}
PARAM_TYPES = {"conifer-file": "inputfile",
               "exomedepth-file": "inputfile",
               "gatk4-file": "inputfile",
               "outdir": "directory",
               "output-prefix": "string"}
TOOL_USAGE = {"conifer_exomedepth": "python comparison.py -t conifer_exomedepth -c conifer_classifcations.txt -e exomedepth_classifications.txt -o comparison_outdir -op conexo_comparison",
              "gatk4_conifer": "python comparison.py -t gatk4_conifer -g gatk4_classifications.txt -c conifer_classifications.txt -o comparison_outdir -op gatcon_comparison",
              "gatk4_exomedepth": "python comparison.py -t gatk4_exomedepth -g gatk4_classifications.txt -e exomedepth_classifications.txt -o comparison_outdir -op gatexo_comparison"}


def main():
    compare_parameters = parpar.get_comparison_parameters(TOOL_CHOICES)
    incorrect_parameters = parpar.parameters_are_ok(compare_parameters, REQUIRED_PARAMS, PARAM_TYPES)

    if len(incorrect_parameters) == 0:
        tool1data = {}
        tool2data = {}
        tool1_label = ""
        tool2_label = ""
        arraydata = ufr.read_array_cnvs(totalsparams["arraycnvs"])

        # Perform comparison between Conifer and ExomeDepth classification data
        if compare_parameters["tool"] == "conifer_exomedepth":
            tool1_label = "Conifer"
            tool2_label = "ExomeDepth"
            print("Implementing...")

        # Perform comparison between GATK4 and Conifer classification data
        if compare_parameters["tool"] == "gatk4_conifer":
            tool1_label = "GATK4"
            tool2_label = "Conifer"
            print("...Reading GATK4 classification data...")
            tool1data = ufr.read_classification_file(compare_parameters["gatk4-file"])
            print("...Reading Conifer classification data...")
            tool2data = ufr.read_classification_file(compare_parameter["conifer-file"])

        # Perform comparison between GATK4 and ExomeDepth classification data
        if compare_parameters["tool"] == "gatk4_exomedepth":
            tool1_label = "GATK4"
            tool2_label = "ExomeDepth"
            print("Implementing...")

        print(f"...Perform the comparison between {tool1_label} and {tool2_label}...")
        comparisondata = comcom.perform_comparison(tool1_label, tool1data, tool2_label, tool2data, arraydata)

        outfilepath = compare_parameters["outdir"] + "/" + compare_parameters["output-prefix"] + ".txt"
        print(f"...Writing comparison data to output file {outfilepath}...")
        wrote_file = ufw.write_comparison_data(outfilepath, comparisondata, tool1_label, tool2_label)
        print(f"...Wrote comparison output file?: {wrote_file}...")
    else:
        print("Please set the following parameters: " + ", ".join(incorrect_parameters))


if __name__ == "__main__":
    main()
    print("DONE!")
