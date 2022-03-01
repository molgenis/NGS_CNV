#!/usr/bin/env python
import os
import sys


def calculate_coverage_ratio(ncfileloc):
    avg_autosomal_coverage = 0
    x_coverage = []
    try:
        with open(ncfileloc, 'r') as ncfile:
            next(ncfile)
            for fileline in ncfile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[0] == "X":
                    if avg_autosomal_coverage == 0:
                        avg_autosomal_coverage = float(filelinedata[5])
                    x_coverage.append(float(filelinedata[4]))
        avg_x_coverage = sum(x_coverage) / len(x_coverage)
        return avg_x_coverage / avg_autosomal_coverage
    except IOError:
        print(f"Could not read {ncfileloc}")


# Collect all the normalized.coverage.txt files
indir = f"{sys.argv[2]}/" if not sys.argv[2].endswith("/") else sys.argv[2]
indirfiles = os.listdir(indir)
coverage_files = [f"{indir}{encf}" for encf in indirfiles if encf.endswith(".normalized.coverage.txt")]


# Make a handy hashmap for the normalized.coverage.txt files
normalized_coverage_files = {}
for ncfile in coverage_files:
    sample_enumber = ncfile.split("/")[-1].split(".")[0]
    normalized_coverage_files[sample_enumber] = ncfile


# Read the sample to sex table
sex_to_bam = []
try:
    with open(sys.argv[1], 'r') as infile:
        next(infile)
        for fileline in infile:
            filelinedata = fileline.strip().split("\t")
            if len(filelinedata) == 3:
                if filelinedata[2] == "UD":
                    sex_to_bam.append(filelinedata[0])
except IOError:
    print("Could not read BAM to sex data")


# Filter the list of normalized.coverage.txt files
for enumber in sex_to_bam:
    if enumber in normalized_coverage_files:
        coverage_ratio = calculate_coverage_ratio(normalized_coverage_files[enumber])
        if coverage_ratio < 0.7:
            print(f"{enumber}: M ({coverage_ratio})")
        elif 0.85 < coverage_ratio < 1.3:
            print(f"{enumber}: F ({coverage_ratio})")
        else:
            print(f"{enumber}: UD ({coverage_ratio})")


# Arguments
# sys.argv[1]: Table linking sample to sex
# sys.argv[2]: Directory with normalized.coverage.txt files
