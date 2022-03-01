#!/usr/bin/env python
import os
import argparse
from ccrscall import CcrsCall


def get_params():
    """Define, receive and return set CLI parameter values."""
    combine_ccrs_args = argparse.ArgumentParser()
    combine_ccrs_args.add_argument("-1", "--nonud-indir", type=str, dest="nonud-indir", required=True, help="Path to directory with combined CCRS ")
    combine_ccrs_args.add_argument("-2", "--ud-indir", type=str, dest="ud-indir", required=True, help="Path to directory with combined CCRS files for UD samples")
    combine_ccrs_args.add_argument("-o", "--output-dir", type=str, dest="output-dir", required=True, help="")
    return vars(combine_ccrs_args.parse_args())


def read_combined_ccrs(infileloc):
    """Read the combined CCRS file."""
    ccrs_calls = {}
    try:
        with open(infileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")

                # Add the sample to the CCRS calls data
                if filelinedata[0] not in ccrs_calls:
                    ccrs_calls[filelinedata[0]] = []
                ccrs_calls[filelinedata[0]].append(CcrsCall(filelinedata[0], filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), int(filelinedata[4]), filelinedata[5], float(filelinedata[6])))
    except IOError:
        print("Could not read combined CCRS file")
    finally:
        return ccrs_calls


def write_ccrs_seg_file(samplenames, batchdata, outfileloc):
    """Write the combined CCRS calls file with non-UD and UD samples, sorted by sample name."""
    try:
        file_written = False
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tChromosome\tStart\tEnd\tNum_Probes\tCall\tSegment_Mean\n")
            for samplename in samplenames:
                for ccrscall in batchdata[samplename]:
                    outfile.write(f"{ccrscall.ccrs_sample}\t{ccrscall.ccrs_chrom}\t{ccrscall.ccrs_start}\t{ccrscall.ccrs_end}\t{ccrscall.ccrs_numofprobes}\t{ccrscall.ccrs_call}\t{ccrscall.ccrs_segmentmean}\n")
        file_written = True
    except IOError:
        print("Could not write non-UD and UD combined CCRS calls file")
    finally:
        return file_written


def main():
    """Combine non-UD and UD sample CCRS calls, sorted by sample name."""
    # Get the parameters.
    combine_ccrs_params = get_params()

    # Make the output directory
    nonud_dir = combine_ccrs_params["nonud-indir"]+"/" if not combine_ccrs_params["nonud-indir"].endswith("/") else combine_ccrs_params["nonud-indir"]
    ud_dir = combine_ccrs_params["ud-indir"]+"/" if not combine_ccrs_params["ud-indir"].endswith("/") else combine_ccrs_params["ud-indir"]
    output_dir = combine_ccrs_params["output-dir"]+"/" if not combine_ccrs_params["output-dir"].endswith("/") else combine_ccrs_params["output-dir"]

    # Collect the combined CCRS call files per batch name
    nonud_ccrs_files = {x.split(".")[0]:x for x in os.listdir(nonud_dir) if x.endswith(".seg")}
    ud_ccrs_files = {x.split(".")[0]:x for x in os.listdir(ud_dir) if x.endswith(".seg")}
    nonud_not_in_ud = nonud_ccrs_files.keys() - ud_ccrs_files.keys()
    ud_not_in_nonud = ud_ccrs_files.keys() - nonud_ccrs_files.keys()

    # Star combining non-UD and UD combined CCRS call files for each batch (BED file).
    print("[-START COMBINING NON-UD AND UD SAMPLES DATA FOR EACH BATCH-]")
    for batchname in nonud_ccrs_files.keys() & ud_ccrs_files.keys():
        print(f"...Processing {batchname}...")
        batch_data = read_combined_ccrs(f"{nonud_dir}{nonud_ccrs_files[batchname]}")
        batch_data.update(read_combined_ccrs(f"{ud_dir}{ud_ccrs_files[batchname]}"))
        sample_names = list(batch_data.keys())
        sample_names.sort()
        wrote_ccrs = write_ccrs_seg_file(sample_names, batch_data, f"{output_dir}{batchname}.called.seg")
        print(f"...Wrote combined CCRS file with non-UD and UD?: {wrote_ccrs}...")

    # Check whether there are batches with only non-UD samples and write the combined CCRS calls files
    if len(nonud_not_in_ud) > 0:
        print("[-START WRITING NON-UD SAMPLES DATA FOR EACH BATCH-]")
        for batchname in nonud_not_in_ud:
            print(f"...Processing {batchname}...")
            batch_data = read_combined_ccrs(f"{nonud_dir}{nonud_ccrs_files[batchname]}")
            sample_names = list(batch_data.keys())
            sample_names.sort()
            wrote_ccrs = write_ccrs_seg_file(sample_names, batch_data, f"{output_dir}{batchname}.called.seg")
            print(f"...Wrote combined CCRS file with non-UD?: {wrote_ccrs}...")

    # Check whether there are batches with only UD samples and write the combined CCRS calls files
    if len(ud_not_in_nonud) > 0:
        print("[-START WRITING UD SAMPLES DATA FOR EACH BATCH-]")
        for batchname in ud_not_in_nonud:
            print(f"...Processing {batchname}...")
            batch_data = read_combined_ccrs(f"{ud_dir}{ud_ccrs_files[batchname]}")
            sample_names = list(batch_data.keys())
            sample_names.sort()
            wrote_ccrs = write_ccrs_seg_file(sample_names, batch_data, f"{output_dir}{batchname}.called.seg")
            print(f"...Wrote combined CCRS file with UD?: {wrote_ccrs}...")


if __name__ == "__main__":
    main()
