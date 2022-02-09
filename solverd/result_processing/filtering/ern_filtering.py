import argparse
from ccrscall import CcrsCall
from erngene import ErnGene
from read_combined_ccrs import read_combined_ccrs


ERN_TO_FILE = {"ERN GENTURIS": "ERN-GENTURIS.2021-06-03.230genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv",
               "ERN ITHACA": "ERN-ITHACA.2021-06-23.3081genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv",
               "ERN-RND": "ERN-RND.2021-07-13.1820genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv",
               "ERN EURO-NMD": "MuscleGeneTable.2021-05-26.611genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv"}
CCRS_HEADER_LINE = ""
global_ccrs_header_line = ""

def get_params():
    """Define, receive and return parameter values."""
    ernfilter_args = argparse.ArgumentParser()
    ernfilter_args.add_argument("-c", "--combine-ern", action="store_true", dest="combine-ern", help="Use the ERN files combined rather than separately to filter")
    ernfilter_args.add_argument("-i", "--infile", type=str, dest="infile", required=True, help="Path to combined CCRS file")
    ernfilter_args.add_argument("-e", "--erndir", type=str, dest="erndir", required=True, help="Path to ERN directory with files to use as filter")
    ernfilter_args.add_argument("-o", "--outdir", type=str, dest="outdir", required=True, help="Path to write output files to")
    ernfilter_args.add_argument("-op", "--out-prefix", type=str, dest="out-prefix", required=True, help="Prefix to use for the output files")
    ernfilter_args.add_argument("-p", "--percentage-overlap", type=str, dest="percentage-overlap", default=25, help="Amount of percentage overlap required")
    ernfilter_args.add_argument("-s", "--samples-to-ern", type=str, dest="samples-to-ern", help="Path to RD3 smapes to ERN file")
    return vars(ernfilter_args.parse_args())


def read_samples_to_ern(samplestoernfileloc):
    """Read the samples to ERN file.

    Parameters
    ----------
    samplestoernfileloc : str
        Path to the samples to ern file
    """
    samples_to_ern = {}
    try:
        with open(samplestoernfileloc, 'r') as samplestoernfile:
            for fileline in samplestoernfile:
                filelinedata = fileline.strip().split("\t")
                samples_to_ern[filelinedata[0]] = filelinedata[1]
    except IOError:
        print("Could not read samples to ERN file")
    finally:
        return samples_to_ern


def read_ern_file(ernfileloc):
    """Read the ERN gene list.

    Parameters
    ----------
    ernfileloc : str
        Path to the ERN genelist

    Returns
    -------
    ern_data : dict
        ERN genes per chromosome
    """
    ern_data = {}
    try:
        with open(ernfileloc, 'r') as ernfile:
            for fileline in ernfile:
                filelinedata = fileline.strip().split("\t")

                # Add the ERN gene to the data
                if filelinedata[1] not in ern_data:
                    ern_data[filelinedata[1]] = []
                ern_data[filelinedata[1]].append(ErnGene(filelinedata[0], filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), int(filelinedata[4]), int(filelinedata[5]), int(filelinedata[6]), int(filelinedata[8])))
    except IOError:
        print("Could not read ERN file")
    finally:
        return ern_data


def determine_ccrs_ern_overlaps(ccrsdata, erndata, samplestoern):
    """Determine which ERN genes overlap with CCRS calls."""
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                if samplename in samplestoern:
                    if samplestoern[samplename] == "UDN-Spain":
                        if ccrscall.ccrs_chrom in erndata["ERN ITHACA"]:
                            ccrscall = get_ern_genes_for_ccrscall(ccrscall, erndata["ERN ITHACA"][ccrscall.ccrs_chrom])
                    elif samplestoern[samplename] in erndata:
                        if ccrscall.ccrs_chrom in erndata[samplestoern[samplename]]:
                            ccrscall = get_ern_genes_for_ccrscall(ccrscall, erndata[samplestoern[samplename]][ccrscall.ccrs_chrom])
    return ccrsdata


def get_ern_genes_for_ccrscall(ccrscall, ernchromdata):
    """Get the overlapping ERN genes for a single CCRS call."""
    for erngene in ernchromdata:
        if ccrs_erngene_overlap(ccrscall, erngene):
            ccrscall.ern_genes.append(erngene)
    if len(ccrscall.ern_genes) > 0:
        ccrscall.keep_call = True
    return ccrscall


def ccrs_erngene_overlap(ccrscall, erngene):
    """Determine whether there is an overlap between the CCRS call and Conrad CNV."""
    if ccrscall.ccrs_chrom == erngene.ern_chrom:
        return ccrscall.ccrs_start <= erngene.ern_padded_stop and erngene.ern_padded_start <= ccrscall.ccrs_end
    return False


def filter_ccrs_with_ern(outfileloc, ccrsdata, ccrsheader, samplenames):
    """Filter the CCRS calls by only keeping the calls overlapping with one or more ERN genes."""
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            #outfile.write(global_ccrs_header_line)
            outfile.write(ccrsheader)
            # outfile.write("Sample\tChromosome\tStart\tEnd\tNum_Probes\tCall\tSegment_Mean\n")
            for samplename in samplenames:
                for chromname in ccrsdata[samplename]:
                    for ccrscall in ccrsdata[samplename][chromname]:
                        if ccrscall.keep_call:
                            outfile.write(ccrscall.to_ccrs_file_line_2())
                            # outfile.write(f"{samplename}\t{ccrscall.ccrs_chrom}\t{ccrscall.ccrs_start}\t{ccrscall.ccrs_end}\t{ccrscall.ccrs_numofprobes}\t{ccrscall.ccrs_call}\t{ccrscall.ccrs_segmentmean}\n")
        file_written = True
    except IOError:
        print("Could not write filtered CCRS segment file")
    finally:
        return file_written


def write_kept_ccrs_calls(outfileloc, ccrsdata, samplenames):
    """Write the CCRS calls that were kept."""
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tCCRS_Call\tCall_Type\tERN_Genes\n")
            for samplename in samplenames:
                for chromname in ccrsdata[samplename]:
                    for ccrscall in ccrsdata[samplename][chromname]:
                        if ccrscall.keep_call:
                            outfile.write(f"{samplename}\t{ccrscall.ccrs_chrom}:{ccrscall.ccrs_start}-{ccrscall.ccrs_end}\t{ccrscall.ccrs_call}\t"+",".join(ccrscall.get_ern_genes())+"\n")
        file_written = True
    except IOError:
        print("Could not write file for kept CCRS calls.")
    finally:
        return file_written


def write_removed_ccrs_calls(outfileloc, ccrsdata, samplenames):
    """Write the CCRS calls that were removed."""
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tCCRS_Call\tCall_Type\n")
            for samplename in samplenames:
                for chromname in ccrsdata[samplename]:
                    for ccrscall in ccrsdata[samplename][chromname]:
                        if not ccrscall.keep_call:
                            outfile.write(f"{samplename}\t{ccrscall.ccrs_chrom}:{ccrscall.ccrs_start}-{ccrscall.ccrs_end}\t{ccrscall.ccrs_call}\n")
        file_written = True
    except IOError:
        print("Could not write file for removed CCRS calls.")
    finally:
        return file_written


def determine_ccrs_combinederd_overlaps(ccrsdata, erndata):
    """Determine which ERN genes of the combined ERN gene lists overlap with CCRS calls.

    Parameters
    ----------
    ccrsdata : dict
    erndata : dict
    """
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                if ccrscall.ccrs_chrom in erndata:
                    ccrscall = get_ern_genes_for_ccrscall(ccrscall, erndata[ccrscall.ccrs_chrom])
    return ccrsdata


def main():
    """Do the main work."""
    ern_filter_params = get_params()
    outdir = ern_filter_params["outdir"]+"/" if not ern_filter_params["outdir"].endswith("/") else ern_filter_params["outdir"]
    outprefix = ern_filter_params["out-prefix"]
    erndir = ern_filter_params["erndir"]+"/" if not ern_filter_params["erndir"].endswith("/") else ern_filter_params["erndir"]

    # Read the CCRS data
    print("[-READING CCRS DATA-]")
    ccrs_data = read_combined_ccrs(ern_filter_params["infile"])
    ccrs_header = ccrs_data[0]
    ccrs_calls = ccrs_data[1]
    print(f"...CCRS DATA: {len(ccrs_data)}...")

    # Read the Samples to ERN table
    print("[-READING SAMPLES TO ERN FILE-]")
    samples_to_ern = read_samples_to_ern(ern_filter_params["samples-to-ern"])
    print(f"...Samples to ERN: {len(samples_to_ern)}...")

    # Read the ERN data
    print("[-READING ERN GENE DATA-]")
    ern_data = {}

    if ern_filter_params["combine-ern"]:
        for ernname in ERN_TO_FILE:
            ern_data.update(read_ern_file(f"{erndir}{ERN_TO_FILE[ernname]}"))
    else:
        for ernname in ERN_TO_FILE:
            ern_data[ernname] = read_ern_file(f"{erndir}{ERN_TO_FILE[ernname]}")
            print(f"...ERN data {ernname}: {len(ern_data[ernname])}...")

    # Determine ERN genes overlapping with CCRS calls
    print("[-DETERMINING ERN GENES OVERLAPPING WITH CCRS CALLS-]")
    if ern_filter_params["combine-ern"]:
        ccrs_data = determine_ccrs_combinederd_overlaps(ccrs_calls, ern_data)
    else:
        ccrs_data = determine_ccrs_ern_overlaps(ccrs_calls, ern_data, samples_to_ern)

    # Sort the list of sample names
    sample_names = list(ccrs_calls.keys())
    sample_names.sort()

    # Determine which calls to retain
    print("[-WRITING FILTERED CCRS SEGMENT FILE-]")
    wrote_segfile = filter_ccrs_with_ern(f"{outdir}{outprefix}.called.seg", ccrs_calls, ccrs_header, sample_names)
    print(f"...Wrote filtered .called.seg file?: {wrote_segfile}...")

    print("[-WRITING CCRS CALL THAT WERE KEPT TO FILE-]")
    wrote_kept_file = write_kept_ccrs_calls(f"{outdir}{outprefix}_kept_calls.txt", ccrs_calls, sample_names)
    print(f"...Wrote kept CCRS calls to file?: {wrote_kept_file}...")

    print("[-WRITING CCRS CALL THAT WERE REMOVED TO FILE-]")
    wrote_removed_file = write_removed_ccrs_calls(f"{outdir}{outprefix}_removed_calls.txt", ccrs_calls, sample_names)
    print(f"...Wrote removed CCRS calls to file?: {wrote_removed_file}...")
    print("[-FINISHED ERN GENE FILTERING-]")


if __name__ == "__main__":
    main()
