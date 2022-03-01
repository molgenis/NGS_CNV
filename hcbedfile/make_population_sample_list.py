#!/usr/bin/env python
import os
import sys

accepted_labels = ["F1", "F2", "M1", "M2"]
label_explanation = {"F1": "Label for population1/female directory",
                     "F2": "Label for population2/female directory",
                     "M1": "Label for population1/male directory",
                     "M2": "Label for population2/male directory"}

if os.path.isdir(sys.argv[1]):
    populationlabel = sys.argv[2].upper()
    if populationlabel in accepted_labels:
        # Collect the BAM files in the directory
        inputdir = f"{sys.argv[1]}/" if not sys.argv[1].endswith("/") else sys.argv[1]
        indirfiles = os.listdir(inputdir)
        bamfiles = [f"{inputdir}{bamfile}" for bamfile in indirfiles if bamfile.endswith(".bam")]

        # Write the population BAM list file
        try:
            with open(sys.argv[3], 'w') as outputfile:
                for bamfile in bamfiles:
                    outputfile.write(f"{bamfile}\t{populationlabel}\n")
        except IOError:
            print(f"Could not write outputfile to {outputfile}")
    else:
        print("Please use one for the following labels:")
        for poplabel in label_explanation:
            print(f"{poplabel} ({label_explanation[poplabel]})")
else:
    print(f"Input directory {sys.argv[1]} does not exist")


# Variables in order:
# sys.argv[1]: Input directory containing BAM files
# sys.argv[2]: Label to use (use F1, F2, M1 or M2)
# sys.argv[3]: Path to write output file to
