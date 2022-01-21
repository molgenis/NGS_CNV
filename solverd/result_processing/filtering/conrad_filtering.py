import os
import argparse
import statistics
from ccrscall import CcrsCall
from conradcnv import ConradCnv
from conradexon import ConradExon
from read_combined_ccrs import read_combined_ccrs


CCRS_TO_CONRAD_CALLTYPE = {'+': "gain", '-': "loss", '0': "na"}

def get_params():
    """Define, receive and return set parameter values."""
    conrad_args = argparse.ArgumentParser()
    conrad_args.add_argument("-a", "--annotate", action="store_true", dest="annotate", help="Annotate CCRS calls with Conrad CNVs")
    conrad_args.add_argument("-m", "--mean-annotation", action="store_true", dest="mean-annotation", help="Annotate CCRS calls with the mean occurrence and frequency")
    conrad_args.add_argument("-i", "--infile", type=str, dest="infile", required=True, help="Path to input directory with CCRS call files")
    conrad_args.add_argument("-o", "--outdir", type=str, dest="outdir", required=True, help="Path to output directory")
    conrad_args.add_argument("-c", "--conrad-cnvs", type=str, dest="conrad-cnvs", required=True, help="Path to the Conrad CNV file")
    conrad_args.add_argument("-e", "--exons-file", type=str, dest="exons-file", required=True, help="Path to file with genenames")
    conrad_args.add_argument("-x", "--xchrom-exons-file", type=str, dest="xchrom-exons-file", help="Path to conrad X chromosome genes")
    conrad_args.add_argument("-p", "--out-prefix", type=str, dest="out-prefix", default="conrad", help="Prefix to use for output files")
    # conrad_args.add_argument("-g", "--genes-file", type=str, dest="genes-file", required=True, help="Path to file with genenames")
    return vars(conrad_args.parse_args())


def read_conrad_cnvs(conradfileloc):
    """Read and return the Conrad CNV data."""
    conradcnvs = {}
    try:
        with open(conradfileloc, 'r') as conradfile:
            next(conradfile)
            for fileline in conradfile:
                filelinedata = fileline.strip().split("\t")

                # Add the Conrad CNV
                if filelinedata[0] not in conradcnvs:
                    conradcnvs[filelinedata[0]] = []
                conradcnvs[filelinedata[0]].append(ConradCnv(filelinedata[0], int(filelinedata[1]), int(filelinedata[2]), int(filelinedata[3]), filelinedata[4], filelinedata[5], filelinedata[6].split("/"), filelinedata[7], float(filelinedata[8])))
    except IOError:
        print(f"Could not read Conrad CNV file {conradfileloc}")
    finally:
        return conradcnvs


def read_conrad_exons(conradexonsfileloc, conradexons):
    """Read the Conrad exon data."""
    try:
        with open(conradexonsfileloc, 'r') as conradexonsfile:
            next(conradexonsfile)
            for fileline in conradexonsfile:
                filelinedata = fileline.strip().split("\t")

                # Add the Conrad exon
                if filelinedata[0] not in conradexons:
                    conradexons[filelinedata[0]] = []
                conradexons[filelinedata[0]].append(ConradExon(filelinedata[0], int(filelinedata[1]), int(filelinedata[2]), filelinedata[3]))
    except IOError:
        print("Could not read Conrad exons file :(")
    finally:
        return conradexons


def determine_ccrs_conradcnv_overlaps(ccrsdata, conraddata):
    """Determine which Conrad CNV's overlap with which CCRS calls."""
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                if chromname in conraddata:
                    ccrscall.conrad_cnvs = get_overlapping_conrad_cnvs(ccrscall, conraddata[chromname])
    return ccrsdata


def get_overlapping_conrad_cnvs(ccrscall, conradchromcnvs):
    """Return Conrad CNVs overlapping with the CCRS calls (RETURNS DATA TO `determine_ccrs_conradcnv_overlaps()`)."""
    conradcnvs = []
    for conradcnv in conradchromcnvs:
        if ccrs_conrad_overlap(ccrscall, conradcnv):
            conradcnvs.append(conradcnv)
    return conradcnvs


def ccrs_conrad_overlap(ccrscall, conradcnv):
    """Determine whether there is an overlap between the CCRS call and Conrad CNV."""
    if ccrscall.ccrs_chrom == conradcnv.cnv_chrom:
        return ccrscall.ccrs_start <= conradcnv.cnv_end and conradcnv.cnv_start <= ccrscall.ccrs_end
    return False


def ccrscall_conradexon_overlap(ccrscall, conradexon):
    """Determine whether there is an overlap between the CCRS call and Conrad exon."""
    if ccrscall.ccrs_chrom == conradexon.exon_chrom:
        return ccrscall.ccrs_start <= conradexon.exon_end and conradexon.exon_start <= ccrscall.ccrs_end
    return False


def conradcnv_conradexon_overlap(conradcnv, conradexon):
    """Determine whether there is an overlap between the Conrad CNV and Conrad exon."""
    if conradcnv.cnv_chrom == conradexon.exon_chrom:
        return conradcnv.cnv_start <= conradexon.exon_end and conradexon.exon_start <= conradcnv.cnv_end
    return False


def add_overlapping_exons_to_ccrs(ccrsdata, exondata):
    """Add overlapping Conrad exons to the CCRS calls that have overlapping Conrad CNVs. Also adds overlapping exons to Conrad CNVs."""
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                if len(ccrscall.conrad_cnvs) > 0:
                    if chromname in exondata:
                        exonlist = get_exons_for_ccrs(ccrscall, exondata[chromname])
                        ccrscall.conrad_exons = exonlist

                        # Add exons to overlapping Conrad CNVs
                        for conradcnv in ccrscall.conrad_cnvs:
                            cexonlist = get_exons_for_conradcnv(conradcnv, exondata[chromname])
                            conradcnv.cnv_exons = cexonlist
    return ccrsdata


def get_exons_for_ccrs(ccrscall, exonchromdata):
    """Get overlapping Conrad exons for a single CCRS call."""
    overlapping_exons = []
    for cexon in exonchromdata:
        if ccrscall_conradexon_overlap(ccrscall, cexon):
            overlapping_exons.append(cexon)
    return overlapping_exons


def get_exons_for_conradcnv(conradcnv, exonchromdata):
    """Get overlapping Coonrad exons for a single Conrad CNV."""
    overlapping_exons = []
    for cexon in exonchromdata:
        if conradcnv_conradexon_overlap(conradcnv, cexon):
            overlapping_exons.append(cexon)
    return overlapping_exons


def filter_ccrs_with_conrad_cnvs(ccrsoutloc, filteredoutloc, ccrsdata):
    """Filter the CCRS data."""
    ordered_sample_names = list(ccrsdata.keys())
    ordered_sample_names.sort()
    wrote_outfiles = False

    try:
        ccrsoutfile = open(ccrsoutloc, 'w')
        filteroutfile = open(filteredoutloc, 'w')

        # Write the headers for both files.
        ccrsoutfile.write("Sample\tChromosome\tStart\tEnd\tNum_Probes\tCall\tSegment_Mean\n")
        filteroutfile.write("Sample\tGATK_Call\tConrad_CNV\tGenes\n")

        # Start filtering the CCRS call using the overlapping Conrad CNVs
        for samplename in ordered_sample_names:
            for chromname in ccrsdata[samplename]:
                for ccrscall in ccrsdata[samplename][chromname]:
                    if len(ccrscall.conrad_cnvs) > 0:
                        # Check that the CCRS call and all overlapping Conrad CNV calls are of the same type and all CCRS genes are in the overlapping Conrad CNVs
                        if ccrscall.can_be_conrad_filtered(CCRS_TO_CONRAD_CALLTYPE):
                            ccrscall.keep_call = False
                            filteroutfile.write(f"{samplename}\t{ccrscall.ccrs_chrom}:{ccrscall.ccrs_start}-{ccrscall.ccrs_end}\t" +",".join(ccrscall.get_conrad_cnv_strings())+ "\t" +",".join(ccrscall.get_gene_names_2())+ "\n")
                        else:
                            ccrscall.keep_call = True
                            ccrsoutfile.write(f"{samplename}\t{ccrscall.ccrs_chrom}\t{ccrscall.ccrs_start}\t{ccrscall.ccrs_end}\t{ccrscall.ccrs_numofprobes}\t{ccrscall.ccrs_call}\t{ccrscall.ccrs_segmentmean}\n")
                    else:
                        ccrscall.keep_call = True
                        ccrsoutfile.write(f"{samplename}\t{ccrscall.ccrs_chrom}\t{ccrscall.ccrs_start}\t{ccrscall.ccrs_end}\t{ccrscall.ccrs_numofprobes}\t{ccrscall.ccrs_call}\t{ccrscall.ccrs_segmentmean}\n")
        ccrsoutfile.close()
        filteroutfile.close()
        wrote_outfiles = True
    except IOError:
        print("Could not filter CCRS call with Conrad CNVs :(")
    finally:
        return wrote_outfiles


def get_num_of_ccrs_with_conradcnvs(ccrsdata):
    """."""
    num_of_ccrs = 0
    total_ccrs_calls = 0
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                total_ccrs_calls += 1
                if len(ccrscall.conrad_cnvs) > 0:
                    num_of_ccrs += 1
    return f"{num_of_ccrs}/{total_ccrs_calls}"


def get_num_of_filtered_ccrscalls(ccrsdata):
    number_of_filtered = 0
    total_ccrs_calls = 0
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                total_ccrs_calls += 1
                if ccrscall.keep_call == False:
                    number_of_filtered += 1
    return f"{number_of_filtered}/{total_ccrs_calls}"


def main():
    """Do the main work."""
    conrad_params = get_params()

    # Make the outdir path
    outdir = conrad_params["outdir"]+"/" if not conrad_params["outdir"].endswith("/") else conrad_params["outdir"]
    outprefix = conrad_params["out-prefix"]

    # Read the combined CCRS file
    print("[-READING COMBINED CCRS FILE-]")
    ccrs_file_data = read_combined_ccrs(conrad_params["infile"])
    ccrs_header = ccrs_file_data[0]
    ccrs_data = ccrs_file_data[1]
    print(f"...CCRS length: {len(ccrs_data)}...")

    # Read the Conrad CNV data and add the names of overlapping genes
    print("[-READING CONRAD DATA-]")
    conrad_cnvs = read_conrad_cnvs(conrad_params["conrad-cnvs"])
    print(f"...Conrad CNV length: {len(conrad_cnvs)}...")
    conrad_exons = {}
    conrad_exons = read_conrad_exons(conrad_params["exons-file"], conrad_exons)
    conrad_exons = read_conrad_exons(conrad_params["xchrom-exons-file"], conrad_exons)
    print(f"...Conrad Exon length: {len(conrad_exons)}...")

    # Determine which Conrad CNV's overlap with which calls.
    print("[-DETERMINING CCRS AND CONRAD CNV OVERLAP-]")
    ccrs_data = determine_ccrs_conradcnv_overlaps(ccrs_data, conrad_cnvs)
    print(f"...Number of samples with overlapping Conrad CNVs: {get_num_of_ccrs_with_conradcnvs(ccrs_data)}...")

    # Determine overlapping Conrad exons for CCRS calls
    print("[-ADDING OVERLAPPING CONRAD EXONS TO CCRS CALL AND ITS OVERLAPPING CONRAD CNVS-]")
    ccrs_data = add_overlapping_exons_to_ccrs(ccrs_data, conrad_exons)

    if conrad_params["annotate"]:
        # Determine the frequency with Conrad CNVs with the lennart method
        print("[-DETERMINING THE OCCURRENCE AND FREQUENCY-]")
        ccrs_data = determine_ccrs_conrad_frequencies_lennart(ccrs_data)

        # Write the annotated calls to file
        print("[-WRIITNG CONRAD ANNOTATED CALLS TO FILE-]")
        wrote_file = write_conrad_annotation_lennart(f"{outdir}{outprefix}.called.seg", ccrs_header, ccrs_data)
    else:
        # Filter the CCRS calls
        print("[-FILTERING CCRS CALLS-]")
        filter_ccrs_with_conrad_cnvs(f"{outdir}{outprefix}.called.igv.seg", f"{outdir}{outprefix}_filtered_ccrs_calls.txt", ccrs_data)
        print(f"...Number of filtered samples: {get_num_of_filtered_ccrscalls(ccrs_data)}...")
        print("[-FINISHED CONRAD FILTERING-]")


def annotate_ccrs_calls_mean(outfileloc, ccrsheader, ccrsdata):
    """Annotate the CCRS calls with mean Conrad CNV occurrences and frequencies."""
    file_written = False
    sample_names = list(ccrsdata.keys())
    sample_names.sort()
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(f"{ccrsheader.strip()}\tconrad_occurrence\tconrad_frequency\n")
            for samplename in sample_names:
                for chromname in ccrsdata[samplename]:
                    for ccrscall in ccrsdata[samplename][chromname]:
                        mean_occurrence = ccrscall.determine_mean_conrad_occurrence()
                        mean_frequency = ccrscall.determine_mean_conrad_frequency()
                        outfile.write(f"{ccrscall.to_ccrs_file_line().strip()}\t{mean_occurrence}\t{mean_frequency}\n")
        file_written = True
    except IOError:
        print("Could not write mean Conrad CNV annotated CCRS file")
    finally:
        return file_written


def annotate_ccrs_calls(outfileloc, ccrsheader, ccrsdata):
    """Annotate the CCRS calls with Conrad CNV occurrences and frequencies."""
    file_written = False
    sample_names = list(ccrsdata.keys())
    sample_names.sort()
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(f"{ccrsheader.strip()}\tconrad_occurrence\tconrad_frequency\n")
            for samplename in sample_names:
                for chromname in ccrsdata[samplename]:
                    for ccrscall in ccrsdata[samplename][chromname]:
                        conrad_occurrences = ccrscall.get_conrad_occurrences()
                        conrad_frequencies = ccrscall.get_conrad_frequencies()
                        outfile.write(f"{ccrscall.to_ccrs_file_line.strip()}\t{conrad_occurrences}\t{conrad_frequencies}\n")
    except IOError:
        print("Could not write Conrad CNV annotated CCRS file")
    finally:
        return file_written


def determine_ccrs_conrad_frequencies_lennart(ccrsdata):
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                c_occur_freq = ccrscall.get_conrad_occurrence_and_frequency(CCRS_TO_CONRAD_CALLTYPE)
                ccrscall.conrad_occurrence = [c_occur_freq[0]]
                ccrscall.conrad_frequency = [c_occur_freq[1]]
    return ccrsdata


def determine_ccrs_conrad_frequencies_lennart_2(ccrsdata):
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromnane]:
                if ccrscall.can_be_conrad_filtered():
                    c_occur_freq = ccrscall.get_max_conrad_occurrence_and_frequency()
                    ccrscall.conrad_occurrence = [c_occurr_freq[0]]
                    ccrscall.conrad_frequency = [c_occurr_freq[1]]
                else:
                    ccrscall.conrad_occurrence = ["0/40"]
                    ccrscall.conrad_frequency = [0.0]
    return ccrsdata


def write_conrad_annotation_lennart(outfileloc, ccrsheader, ccrsdata):
    """."""
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(f"{ccrsheader.strip()}\tConrad_Occurrence\tConrad_Frequency\n")
            for samplename in ccrsdata:
                for chromname in ccrsdata[samplename]:
                    for ccrscall in ccrsdata[samplename][chromname]:
                        outfile.write(ccrscall.to_ccrs_file_line())
        file_written = True
    except IOError:
        print("Could not write Conrad annotated file :|")
    finally:
        return file_written


if __name__ == "__main__":
    main()
