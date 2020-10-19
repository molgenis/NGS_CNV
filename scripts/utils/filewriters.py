def write_output_file(outfilepath, headerline, intervaldata):
    """Write intervals overlapping with a genomic region to an output file.

    Parameters
    ----------
    outfilepath : str
        Path to the output file.
    intervaldata : dict
        Overlapping intervals to write.
    """
    try:
        with open(outfilepath, "w") as outfile:
            outfile.write(headerline)
            for region in intervaldata:
                for overlapinterval in intervaldata[region]:
                    interval_genomicpos = overlapinterval.split("\t")[0:3]
                    intervallabel = make_interval_label(interval_genomicpos[0], interval_genomicpos[1], interval_genomicpos[2])
                    outfile.write(f"{overlapinterval}\t{region}\t{intervallabel}\n")
    except IOError:
        print(f"Could not write to output file {outfilepath}")


def write_seg_output_file(outfilepath, intervaldata):
    """Writes the output file for a SEG input file.
    
    Parameters
    ----------
    outfilepath : str
        Path to write output file to.
    intervaldata : dict
        Intervals to write to output file.
    """
    header_to_add = "CONTIG\tSTART\tEND\tNUM_POINTS_COPY_RATIO\tMEAN_LOG2_COPY_RATIO\tCALL\tOVERLAP\tLABEL\n"
    write_output_file(outfilepath, header_to_add, intervaldata)


def write_tsv_output_file(outfilepath, intervaldata):
    """Writes the output file for a TSV input file.
    
    Parameters
    ----------
    outfilepath : str
        Path to write output file to.
    intervaldata : dict
        Intervals to write to output file.
    """
    header_to_add = "CONTIG\tSTART\tEND\tLOG2_COPY_RATIO\tOVERLAP\tLABEL\n"
    write_output_file(outfilepath, header_to_add, intervaldata)


def write_cac_output_file(outfilepath, intervaldata):
    """Write intervals overlapping with a genomic region to an output file.

    Parameters
    ----------
    outfilepath : str
        Path to the output file.
    intervaldata : dict
        Overlapping intervals to write.
    """
    try:
        with open(outfilepath, "w") as outfile:
            outfile.write("CONTIG\tPOSITION\tREF_COUNT\tALT_COUNT\tREF_NUCLEOTIDE\tALT_NUCLEOTIDE\tOVERLAP\tLABEL\n")
            for region in intervaldata:
                for overlapinterval in intervaldata[region]:
                    interval_genomicpos = overlapinterval.split("\t")[0:2]
                    intervallabel = make_interval_label(interval_genomicpos[0], interval_genomicpos[1], interval_genomicpos[1])
                    outfile.write(f"{overlapinterval}\t{region}\t{intervallabel}\n")
    except IOError:
        print(f"Could not write to output file {outfilepath}")


def write_missedfound_arraycnvs(missedfound_arraycnvs, outfileloc):
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tArray_CNV\t#_Probes\n")
            for samplename in missedfound_arraycnvs:
                for mfa_cnv in missedfound_arraycnvs[samplename]:
                    outfile.write(f"{samplename}\t{mfa_cnv.get_region()}\t{mfa_cnv.num_of_probes()}\n")
        file_written = True
    except IOError:
        print(f"Could not write missed/found array cnv data to output file: {outfileloc}")
    finally:
        return file_written


def write_missedfound_summary(missedfound_summary, outfileloc):
    """Write the summary data for missed/found array CNVs.

    Parameters
    ----------
    missedfound_summary : dict
        Missed/Found summary data to write
    outfileloc: str
        Path to write output file to
    """
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Cnv_Type\tCount\n")
            for cnvtype in missedfound_summary:
                outfile.write(f"{cnvtype}\t{missedfound_summary[cnvtype]}\n")
        file_written = True
    except IOError:
        print(f"Could not write missed/found summary results to output file: {outfileloc}")
    finally:
        return file_written


def write_classification_totals(totalsdata, outfileloc):
    wrote_file = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Classification\tTotal\n")
            for totalslabel in totalsdata:
                outfile.write(f"{totalslabel}\t{totalsdata[totalslabel]}\n")
        wrote_file = True
    except IOError:
        print(f"Could not write classification totals to {outfileloc}")
    finally:
        return wrote_file


def write_duplicate_fpregions(outfileloc, dupfpregions):
    """Write duplicate False Positive regions to an output file.

    Parameters
    ----------
    outfileloc : str
        Path to write duplicate False Positive regions to
    dupfpregions : dict
        False Positive regions and counts
    """
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Region\tCount\n")
            for dupregion in dupfpregions:
                outfile.write(f"{dupregion}\t{dupfpregions[dupregion]}\n")
    except IOError:
        print(f"Could not write duplicate False Positive regions to {outfileloc}")


def write_unique_fpregions(outfileloc, unifpregions):
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Region\tCount\n")
            for uniregion in unifpregions:
                outfile.write(f"{uniregion}\t{unifpregions[uniregion]}\n")
    except IOError:
        print(f"Could not write unique False Positive regions to {outfileloc}")


def write_similar_fpregions(outfileloc, simfpregions):
    """Write similar False Positive regions to an output file.

    Parameters
    ----------
    outfileloc : str
        Path to write similar False Positive regions to
    simfpegions : dict
        Similar False Positive regions and overlapping regions
    """
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Region\tOverlap_Region\tOverlap_Percentage\n")
            for simregion in simfpregions:
                for overlapregion in simfpregions[simregion].overlaps:
                    outfile.write(f"{simregion}\t{overlapregion}\t{simfpregions[simregion].overlaps[overlapregion]}\n")
    except IOError:
        print(f"Could not write similar False Positive region to {outfileloc}")
