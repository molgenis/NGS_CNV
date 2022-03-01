#!/usr/bin/env python
import argparse


def get_params():
    """Define, receive and return set parameter values."""
    ccrs_probes_args = argparse.ArgumentParser()
    ccrs_probes_args.add_argument("-i", "--infile", type=str, dest="infile", required=True, help="")
    ccrs_probes_args.add_argument("-o", "--outdir", type=str, dest="outdir", required=True, help="")
    ccrs_probes_args.add_argument("-p", "--prefix", type=str, dest="prefix", required=True, help="")
    return vars(ccrs_probes_args.parse_args())


def read_ccrs_probes(ccrsfileloc):
    """Read the CCRS data and return the numbers of probes per sample."""
    call_probes = {}
    try:
        with open(ccrsfileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")

                # Add the call type (dup or del) to the dictionary
                if filelinedata[5] not in call_probes:
                    call_probes[filelinedata[5]] = {}

                # Add the samplename to the dictionary
                if filelinedata[0] not in call_probes[filelinedata[5]]:
                    call_probes[filelinedata[5]][filelinedata[0]] = []

                # Add the actual number of probes
                call_probes[filelinedata[5]][filelinedata[0]].append(filelinedata[4])
    except IOError:
        print("Could not read combined CCRS file")
    finally:
        return call_probes


def write_probes_per_call(outfileloc, ccrsprobes):
    """Write the number of probes per sample to file."""
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tProbes\n")
            for samplename in ccrsprobes:
                for probenum in ccrsprobes[samplename]:
                    outfile.write(f"{samplename}\t{probenum}\n")
        file_written = True
    except IOError:
        print("Could not write probes file")
    finally:
        return file_written


def main():
    """Do the main work."""
    ccrs_probes_params = get_params()
    outdir = ccrs_probes_params["outdir"]+"/" if not ccrs_probes_params["outdir"].endswith("/") else ccrs_probes_params["outdir"]
    dup_out_loc = f"{outdir}" +ccrs_probes_params["prefix"]+ "_dup_probes.txt"
    del_out_loc = f"{outdir}" +ccrs_probes_params["prefix"]+ "_del_probes.txt"

    # Process the data
    ccrs_probe_data = read_ccrs_probes(ccrs_probes_params["infile"])
    wrote_dup_probes = write_probes_per_call(dup_out_loc, ccrs_probe_data['+'])
    wrote_del_probes = write_probes_per_call(del_out_loc, ccrs_probe_data['-'])

    # Check whether the output file have been written
    print(f"Wrote duplication probes output file?: {wrote_dup_probes}")
    print(f"Wrote deletion probes output file?: {wrote_del_probes}")


if __name__ == "__main__":
    main()
