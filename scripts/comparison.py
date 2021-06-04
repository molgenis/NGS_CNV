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
TOOL_CHOICES = ["arraycnvs", "false_positives", "true_positives"]
REQUIRED_PARAMS = {"arraycnvs": ["arraycnvs", "file1", "file2", "label1", "label2", "outdir", "output-prefix"],
                   "false_positives": ["file1", "file2", "label1", "label2", "outdir", "output-prefix"],
                   "true_positives": ["file1", "file2", "label1", "label2", "outdir", "output-prefix"]}
OPTIONAL_PARAMS = {}
PARAM_TYPES = {"arraycnvs": "inputfile",
               "file1": "inputfile",
               "file2": "inputfile",
               "label1": "string",
               "label2": "string",
               "outdir": "directory",
               "output-prefix": "string"}
TOOL_USAGE = {"conifer_exomedepth": "python comparison.py -t conifer_exomedepth -c conifer_classifcations.txt -e exomedepth_classifications.txt -o comparison_outdir -op conexo_comparison",
              "gatk4_conifer": "python comparison.py -t gatk4_conifer -g gatk4_classifications.txt -c conifer_classifications.txt -o comparison_outdir -op gatcon_comparison",
              "gatk4_exomedepth": "python comparison.py -t gatk4_exomedepth -g gatk4_classifications.txt -e exomedepth_classifications.txt -o comparison_outdir -op gatexo_comparison",
              "gatk4_gatk4": "python comparison.py -t gatk4_gatk4 -g gatk4_classifications.txt -g2 gatk4_classifications_2.txt -a array_goldstandard.txt -o comparison_outdir -op gatkgatk_comparison",
              "conifer_conifer": "python comparison.py -t conifer_conifer -c conifer_classifcations.txt -c2 conifer_classifcations_2.txt -a array_goldstandard.txt -o comparison_outdir -op concon_comparison",
              "exomedepth_exomedepth": "python comparison.py -t exomedepth_exomedepth -e exomedepth_classifications.txt -e2 exomedepth_classifications_2.txt -a array_goldstandard.txt -o comparison_outdir -op exdexd_comparison"}


def main():
    compare_parameters = parpar.get_comparison_parameters(TOOL_CHOICES)
    incorrect_parameters = parpar.parameters_are_ok(compare_parameters, REQUIRED_PARAMS, PARAM_TYPES)

    if len(incorrect_parameters) == 0:
        tool1_label = compare_parameters["label1"]
        tool2_label = compare_parameters["label2"]
        print(f"...Reading {tool1_label} classification data...")
        tool1data = ufr.read_classification_file(compare_parameters["file1"])
        print(f"...Reading {tool2_label} classification data...")
        tool2data = ufr.read_classification_file(compare_parameters["file2"])

        # Perform comparison between two tools for found array CNVs
        if compare_parameters["tool"] == "arraycnvs":
            print("...Reading array CNV data...")
            arraydata = ufr.read_array_cnvs(compare_parameters["arraycnvs"])
            print(f"...Perform the comparison between {tool1_label} and {tool2_label}...")
            comparisondata = comcom.perform_comparison(tool1_label, tool1data, tool2_label, tool2data, arraydata, compare_parameters["tp-per-acnv"])

            outfilepath = compare_parameters["outdir"] + "/" + compare_parameters["output-prefix"] + ".txt"
            print(f"...Writing comparison data to output file {outfilepath}...")
            wrote_file = ufw.write_comparison_data(outfilepath, comparisondata, tool1_label, tool2_label)
            print(f"...Wrote comparison output file?: {wrote_file}...")

        # Perform comparison between two tools for False Positives
        if compare_parameters["tool"] == "false_positives":
            comparisondata = comcom.compare_fps(tool1_label, tool1data, tool2_label, tool2data)
            outfilepath = compare_parameters["outdir"] + "/" + compare_parameters["output-prefix"] + ".txt"
            wrote_file = ufw.write_fp_comparison(outfilepath, comparisondata, tool1_label, tool2_label)
            print(f"...Wrote comparison output file?: {wrote_file}...")

        # Perform comparison between two tools for True Positives
        if compare_parameters["tool"] == "true_positives":
            comparisondata = comcom.compare_tps(tool1_label, tool1data, tool2_label, tool2data)
            outfilepath = compare_parameters["outdir"] + "/" + compare_parameters["output-prefix"] + ".txt"
            wrote_file = ufw.write_tp_comparison(outfilepath, comparisondata, tool1_label, tool2_label)
            print(f"...Wrote comparison output file?: {wrote_file}...")
    else:
        print("Please set the following parameters: " + ", ".join(incorrect_parameters))


if __name__ == "__main__":
    main()
    print("DONE!")
