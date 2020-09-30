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
