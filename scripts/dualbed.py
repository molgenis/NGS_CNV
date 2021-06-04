import os
import argparse

# Import required classes
import parameters.parameters as parpar


# Make parameter defining variables
TOOL_CHOICES = ["classification", "combine", "filter"]
REQUIRED_PARAMS = {"classification": ["infile1", "infile2", "outdir", "percent-overlap"],
                   "combine": ["infile1", "infile2", "outfile"],
                   "filter": ["infile1", "outfile"]}
OPTIONAL_PARAMS = {"classification": ["percent-overlap"],
                   "combine": ["no-unique"]}
PARAM_TYPES = {"infile1": "inputfile",
               "infile2": "inputfile",
               "outfile": "outputfile",
               "outdir": "directory",
               "percent-overlap": "integer"}


def main():
    cbd_params = parpar.get_dualbed_parameters(TOOL_CHOICES)
    incorrect_parameters = parpar.parameters_are_ok(cbd_params, REQUIRED_PARAMS, PARAM_TYPES)

    if len(incorrect_parameters) == 0:
        # Perform the dualBED classification
        if cbd_params["tool"] == "classification":
            print("...Implementing...")

        # Combine two dualBED result files
        if cbd_params["tool"] == "combine":
            print("...Start combining two dualBED files...")
            wrote_file = combine_dualbed_files(cbd_params["infile1"], cbd_params["infile2"], cbd_params["outfile"], cbd_params["no-unique"])
            print(f"...Wrote combined output file?: {wrote_file}...")

        # Filter dualBED file
        if cbd_params["tool"] == "filter":
            print("...Implementing...")
    else:
        print("Please set the following parameters: " + ", ".join(incorrect_parameters))


def combine_dualbed_files(dbfile1, dbfile2, outfileloc, nounique):
    """Combine two dualBED input files and write it to an output file

    Parameters
    ----------
    dbfile1 : str
        Path to first dualBED input file
    dbfile2 : str
        Path to second dualBED input file
    outfileloc : str
        Path to write output file to
    nounique : bool
        Exclude Unique calls if True, Include if False

    Returns
    -------
    file_written : bool
        True if output file has been succesfully written, False if not
    """
    file_written = False
    try:
        outfile = open(outfileloc, 'w')

        # Write first input file entries to the output file
        print("...Adding first dua;BED file...")
        with open(dbfile1, 'r') as infile1:
            for fileline in infile1:
                add_line = True
                filelinedata = fileline.strip().split("\t")
                if filelinedata[-1] == "Unique" and nounique:
                    add_line = False
                if add_line:
                    outfile.write(fileline)

        # Write second input file entries to the output file
        print("...Adding second dualBED file...")
        with open(dbfile2, 'r') as infile2:
            next(infile2)
            for fileline in infile2:
                add_line = True
                filelinedata = fileline.strip().split("\t")
                if filelinedata[-1] == "Unique" and nounique:
                    add_line = False
                if add_line:
                    outfile.write(fileline)

        outfile.close()
        file_written = True
    except IOError:
        print("Could not combine two dualBED input files")
    finally:
        return file_written


if __name__ == "__main__":
    main()
