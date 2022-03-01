#!/usr/bin/env python
import os
import argparse
from ccrscall import CcrsCall
from read_combined_ccrs import read_combined_ccrs


ERN_TO_FILE = {"ERN GENTURIS": "ERN-GENTURIS.2021-06-03.230genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv",
               "ERN ITHACA": "ERN-ITHACA.2021-06-23.3081genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv",
               "ERN-RND": "ERN-RND.2021-07-13.1820genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv",
               "ERN EURO-NMD": "MuscleGeneTable.2021-05-26.611genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv"}


def get_params():
    """Define, receive and return set CLI parameter values."""
    erncombine_args = argparse.ArgumentParser()
    erncombine_args.add_argument("-c", "--ccrs-dir", type=str, required=True, dest="ccrs-dir", help="Path to directory with combined CCRS files")
    erncombine_args.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", help="Path to output directory")
    erncombine_args.add_argument("-e", "--ern-dir", type=str, required=True, dest="ern-dir", help="Path to directory with ERN gene lists")
    erncombine_args.add_argument("-s", "--samples-to-ern", type=str, required=True, dest="samples-to-ern", help="Path to file containing samples to ERN")
    return vars(erncombine_args.parse_args())


def get_ccrs_calls(ccrsdir, ccrsfiles):
    """Read and return CCRS calls from all files.

    Parameters
    ----------
    ccrsfiles : dict
        Paths to combined CCRS files
    """
    all_ccrs_calls = {}
    batchnames = list(ccrsfiles.keys())
    batchnames.sort()
    ccrsheader = ""

    # Read each combined CCRS file
    for batchname in batchnames:
        print(f"...Reading CCRS call for {batchname}...")
        batchdata = read_combined_ccrs(f"{ccrsdir}{ccrsfiles[batchname]}")
        all_ccrs_calls.update(batchdata[1])
        if ccrsheader == "":
            ccrsheader = batchdata[0]
    return [ccrsheader, all_ccrs_calls]


def read_samples_to_ern(infileloc):
    """Read the sample names for each ERN.

    Parameters
    ----------
    infileloc : str
        Path to samples to ERN file

    Returns
    -------
    ern_samples : dict
        Sample names per ERN
    """
    ern_samples = {}
    try:
        with open(infileloc, 'r') as infile:
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[1] not in ern_samples:
                    ern_samples[filelinedata[1]] = []
                ern_samples[filelinedata[1]].append(filelinedata[0])
    except IOError:
        print("Could not read samples to ERN file")
    finally:
        return ern_samples


def get_calls_per_ern(ernsamples, ccrscalls):
    """Get the actual calls per ERN.

    Parameters
    ----------
    ernsamples : dict
        Sample names per ERN
    ccrscalls : dict
        CCRS calls per sample, per chrom
    """
    ern_calls = {}
    for ernname in ernsamples:
        ern_calls[ernname] = {}
        for samplename in ernsamples[ernname]:
            if samplename in ccrscalls:
                ern_calls[ernname][samplename] = ccrscalls[samplename]
    return ern_calls


def write_ern_ccrs_calls(outdir, ccrsheader, erncalls):
    """Write the CCRS calls per ERN to file.

    Parameters
    ----------
    outdir : str
        Path to output directory
    ccrsheader : str
        Header of the read combined CCRS file
    erncalls : dict
        CCRS calls per ERN
    """
    for ernname in erncalls:
        wrote_file = write_erncalls_file(f"{outdir}{ernname}_calls.txt", ccrsheader, erncalls[ernname])
        print(f"...Wrote CCRS calls file for {ernname}?: {wrote_file}")


def write_erncalls_file(outfileloc, ccrsheader, ernccrscalls):
    """Write CCRS calls for a single ERN.

    Parameters
    ----------
    outfileloc : str
        Path to write output file to
    ccrsheader : str
        Header of the read combined CCRS file
    ernccrscalls : dict
        CCRS calls a single ERN saved per sample, per chrom
    """
    samplenames = list(ernccrscalls.keys())
    samplenames.sort()
    file_written = False

    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(ccrsheader)
            for samplename in samplenames:
                for chromname in ernccrscalls[samplename]:
                    for ccrscall in ernccrscalls[samplename][chromname]:
                        outfile.write(ccrscall.to_ccrs_file_line_2())
        file_written = True
    except IOError:
        print("Could not write ERN CCRS calls file")
    finally:
        return file_written


def main():
    """Do the main work."""
    erncombine_params = get_params()
    # ccrs_files = [x for x in os.listdir(erncombine_params["ccrs-dir"]) if x.endswith(".called.seg")]

    print("[-READING ALL CCRS CALLS-]")
    ccrs_files = {x.split('.')[0]: x for x in os.listdir(erncombine_params["ccrs-dir"]) if x.endswith(".called.seg")}
    ccrs_data = get_ccrs_calls(erncombine_params["ccrs-dir"], ccrs_files)
    ccrs_header = ccrs_data[0]
    ccrs_calls = ccrs_data[1]

    print("[-READING SAMPLES TO ERN-]")
    ern_samples = read_samples_to_ern(erncombine_params["samples-to-ern"])

    print("[-LINKING CCRS CALLS TO ERN-]")
    ern_calls = get_calls_per_ern(ern_samples, ccrs_calls)

    print("[-WRITING CCRS CALLS PER ERN-]")
    write_ern_ccrs_calls(erncombine_params["outdir"], ccrs_header, ern_calls)


if __name__ == "__main__":
    main()
