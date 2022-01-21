import argparse
from ccrscall import CcrsCall
from erngene import ErnGene

ERN_TO_FILE = {"ERN GENTURIS": "ERN-GENTURIS.2021-06-03.230genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv",
               "ERN ITHACA": "ERN-ITHACA.2021-06-23.3081genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv",
               "ERN-RND": "ERN-RND.2021-07-13.1820genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv",
               "ERN EURO-NMD": "MuscleGeneTable.2021-05-26.611genes.HGNC.CHR.Start.End.Len.ModStart.ModEnd.ModLen.tsv"}


def get_params():
    """Define, receive and return parameter values."""
    ernfilter_args = argparse.ArgumentParser()
    ernfilter_args.add_argument("-c", "--combine-ern", action="store_true", dest="combine-ern", help="Use the ERN files combined rather than separately to filter")
    ernfilter_args.add_argument("-i", "--infile", type=str, dest="infile", required=True, help="Path to combined CCRS file")
    ernfilter_args.add_argument("-e", "--erndir", type=str, dest="erndir", required=True, help="Path to ERN directory with files to use as filter")
    ernfilter_args.add_argument("-o", "--outdir", type=str, dest="outdir", required=True, help="Path to write output files to")
    ernfilter_args.add_argument("-op", "--out-prefix", type=str, dest="out-prefix", required=True, help="Prefix to use for the output files")
    ernfilter_args.add_argument("-p", "--percentage-overlap", type=str, dest="percentage-overlap", default=25, help="Amount of percentage overlap required")
    ernfilter_args.add_argument("-s", "--samples-to-ern", type=str, dest="samples-to-ern", help="Path to samples to ERN file")
    return vars(ernfilter_args.parse_args())


def read_combined_ccrs(ccrsfileloc):
    """Read the combined CCRS data."""
    ccrsdata = {}
    try:
        with open(ccrsfileloc, 'r') as ccrsfile:
            global_ccrs_header_line = next(ccrsfile)
            for fileline in ccrsfile:
                filelinedata = fileline.strip().split("\t")

                # Add the sample name and chromosome as dictionary entries
                if filelinedata[0] not in ccrsdata:
                    ccrsdata[filelinedata[0]] = {}
                if filelinedata[1] not in ccrsdata[filelinedata[0]]:
                    ccrsdata[filelinedata[0]][filelinedata[1]] = []
                ccrs_call = CcrsCall(filelinedata[0], filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), int(filelinedata[4]), filelinedata[5], float(filelinedata[6]))

                # Check if the there is frequency annotation
                if len(filelinedata) >= 9:
                    ccrs_call.ccrs_occurrence = filelinedata[7]
                    ccrs_call.ccrs_frequency = float(filelinedata[8])

                # Check if there is conrad annotation
                if len(filelinedata) >= 11:
                    ccrs_call.conrad_occurrence = filelinedata[9]
                    ccrs_call.conrad_frequency = float(filelinedata[10])

                # Check if there is gnomad annotation
                if len(filelinedata) >= 12:
                    ccrs_call.gnomad_frequency = float(filelinedata[11])

                # Add the CCRS call to the CCRS data
                ccrsdata[filelinedata[0]][filelinedata[1]].append(ccrs_call)
    except IOError:
        print("Could not read combined CCRS file :(")
    finally:
        return ccrsdata


def read_samples_to_ern(samplestoernfileloc):
    """Read the samples to ERN file."""
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
    """Read the ERN gene list."""
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


def determine_ccrs_combinederd_overlaps(ccrsdata, erndata):
    """Determine which ERN genes of the combined ERN gene lists overlap with CCRS calls."""
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                if ccrscall.ccrs_chrom in erndata:
                    ccrscall = get_ern_genes_for_ccrscall(ccrscall, erndata[ccrscall.ccrs_chrom])
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


def determine_call_frequencies(ccrsdata):
    """Determine the counts for each CCRS call based on the genes."""
    numofsamples = len(ccrsdata)
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                ccrscallcount = get_ccrs_call_count(ccrscall, samplename, ccrsdata)
                ccrscall.ccrs_occurrence = f"{ccrscallcount}/{numofsamples}"
                ccrscall.ccrs_frequency = round((ccrscallcount/numofsamples)*100, 2)
    return ccrsdata


def get_ccrs_call_count(ccrscall, ccrssample, ccrsdata):
    """Determine the occurrence of each CCRS call."""
    callcount = 1
    for samplename in ccrsdata:
        if samplename != ccrssample:
            for chromname in ccrsdata[samplename]:
                for otherccrscall in ccrsdata[samplename][chromname]:
                    if ccrscall.ccrs_call == otherccrscall.ccrs_call:
                        # Check if all the genes are in the other call
                        if len(set(ccrscall.get_ern_genes()) - set(otherccrscall.get_ern_genes())) == 0:
                            callcount += 1
    return callcount


def ccrs_calls_overlap(ccrscall1, ccrscall2):
    """Determine whether there is an overlap between the CCRS call and Conrad CNV."""
    if ccrscall.ccrs_chrom == erngene.ern_chrom:
        return ccrscall1.ccrs_start <= ccrscall2.ccrs_end and ccrscall2.ccrs_start <= ccrscall1.ccrs_end
    return False


def annotate_ccrs_file(outfileloc, ccrsdata, samplenames):
    """Write the combined CCRS calls file with annotated occurrence and frequency."""
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tChromosome\tStart\tEnd\tNum_Probes\tCall\tSegment_Mean\tCall_Occurrence\tCall_Frequency\n")
            for samplename in samplenames:
                for chromname in ccrsdata[samplename]:
                    for ccrscall in ccrsdata[samplename][chromname]:
                        outfile.write(ccrscall.to_ccrs_file_line())
        file_written = True
    except IOError:
        print("Could not write frequency annotated combined CCRS file")
    finally:
        return file_written


def main():
    """Annotate CCRS calls with frequency."""
    freq_filter_params = get_params()
    outdir = freq_filter_params["outdir"]+"/" if not freq_filter_params["outdir"].endswith("/") else freq_filter_params["outdir"]
    outprefix = freq_filter_params["out-prefix"]
    erndir = freq_filter_params["erndir"]+"/" if not freq_filter_params["erndir"].endswith("/") else freq_filter_params["erndir"]

    # Read the combined CCRS file
    print("[-READING COMBINED CCRS DATA-]")
    ccrs_data = read_combined_ccrs(freq_filter_params["infile"])

    print("[-READING SAMPLE TO ERN TABLE-]")
    samples_to_ern = read_samples_to_ern(freq_filter_params["samples-to-ern"])
    sample_names = list(ccrs_data.keys())
    sample_names.sort()

    # Read the ERN files into a single file
    print("[-READING ERN FILES-]")
    ern_data = {}
    for ernname in ERN_TO_FILE:
        print(f"...Reading ERN file {ernname}...")
        ern_data.update(read_ern_file(f"{erndir}{ERN_TO_FILE[ernname]}"))

    # Determine the ERN gene overlaps
    print("[-DETERMINE CCRS & ERN GENE OVERLAPS-]")
    ccrs_data = determine_ccrs_combinederd_overlaps(ccrs_data, ern_data)

    # Determine the call frequencies
    print("[-DETERMINE THE CCRS CALL OCCURRENCES AND FREQUENCIES-]")
    ccrs_data = determine_call_frequencies(ccrs_data)

    # Write the annotated combined CCRS file
    print("[-WRITING OCCURRENCE AND FREQUENCY ANNOTATED CCRS FILE-]")
    wrote_file = annotate_ccrs_file(f"{outdir}{outprefix}.called.seg", ccrs_data, sample_names)
    print(f"...Wrote frequency annotated CCRS calls file?: {wrote_file}...")


if __name__ == "__main__":
    main()
