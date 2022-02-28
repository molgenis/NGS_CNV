#!/usr/bin/env python
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


def write_comparison_data(outfilepath, comparisondata, tool1_label, tool2_label):
    """Write comparison results to a specified output file.

    Parameters
    ----------
    outfilepath :str
        Path to write output file to
    comparisondata : dict
        Comparison data to write to file

    Returns
    -------
    file_written : bool
        True if file has been written, False if not
    """
    file_written = False
    comparison_points = ["Total calls", "Call types", "Found array CNVs", "Shared array CNVs", "Shared array types", "Unique array CNVs", "Unique array types", "Shared overlap"]

    try:
        with open(outfilepath, "w") as outfile:
            outfile.write(f"Results for comparison between {tool1_label} and {tool2_label}\n\n")

            for comparepoint in comparison_points:
                # Write data for the total call comparison.
                if comparepoint == "Total calls":
                    outfile.write("[-Total calls made-]\n")
                    outfile.write(f"{tool1_label}: {comparisondata[comparepoint][0]}\n")
                    outfile.write(f"{tool2_label}: {comparisondata[comparepoint][1]}\n")
                    outfile.write("\n\n")

                # Write data for the call type comparison
                if comparepoint == "Call types":
                    outfile.write("[-Call types made-]\n")
                    outfile.write(f"{tool1_label} | {tool2_label}\n")
                    for classificationlabel in comparisondata[comparepoint]:
                        outfile.write(f"* {classificationlabel}: {comparisondata[comparepoint][classificationlabel][0]} | {comparisondata[comparepoint][classificationlabel][1]}\n")
                    outfile.write("\n\n")

                # Write data for the shared found array CNVs
                if comparepoint == "Shared array CNVs":
                    num_of_shared_acnvs = 0
                    outfile.write("[-Array CNVs found by both-]\n")
                    for samplename in comparisondata[comparepoint]:
                        if len(comparisondata[comparepoint][samplename]) > 0:
                            num_of_shared_acnvs += len(comparisondata[comparepoint][samplename])
                            outfile.write("* [" + samplename + "]: " + ", ".join(comparisondata[comparepoint][samplename]) + "\n")
                    outfile.write(f"{tool1_label} and {tool2_label} both identified {num_of_shared_acnvs} array CNVs\n")
                    outfile.write("\n\n")

                # Write data for thje shared array CNV types
                if comparepoint == "Shared array types":
                    outfile.write("Array CNV types found by both:\n")
                    for acnvtype in comparisondata[comparepoint]:
                        outfile.write(f"{acnvtype}: {comparisondata[comparepoint][acnvtype]}\n")
                    outfile.write("\n\n")

                # Write data for unique found array CNVs
                if comparepoint == "Unique array CNVs":
                    outfile.write("[-Unique array CNVs found per tool-]\n")

                    # Write the unique array CNVs for tool 1
                    outfile.write(f"Unique array CNVs found by {tool1_label}\n")
                    num_of_tool1_acnvs = 0
                    for samplename in comparisondata[comparepoint]["tool1"]:
                        num_of_tool1_acnvs += len(comparisondata[comparepoint]["tool1"][samplename])
                        outfile.write("* " + samplename + ": " + ", ".join(comparisondata[comparepoint]["tool1"][samplename]) + "\n")
                    outfile.write(f"{tool1_label} found {num_of_tool1_acnvs} unique array CNVs\n\n")

                    # Write the unique array CNVs for tool 2
                    outfile.write(f"Unique array CNVs found by {tool2_label}\n")
                    num_of_tool2_acnvs = 0
                    for samplename in comparisondata[comparepoint]["tool2"]:
                        num_of_tool2_acnvs += len(comparisondata[comparepoint]["tool2"][samplename])
                        outfile.write("* " + samplename + ": " + ", ".join(comparisondata[comparepoint]["tool2"][samplename]) + "\n")
                    outfile.write(f"{tool2_label} found {num_of_tool2_acnvs} unique array CNVs\n")
                    outfile.write("\n\n")

                # Write data for unique array CNV types found per tool
                if comparepoint == "Unique array types":
                    outfile.write("[-Unique array CNV types per tool-]\n")
                    outfile.write(f"Unique array types for {tool1_label}\n")
                    for acnvtype in comparisondata[comparepoint]["tool1"]:
                        outfile.write(acnvtype + ": " + str(comparisondata[comparepoint]["tool1"][acnvtype]) + "\n")
                    outfile.write("\n")

                    outfile.write(f"Unique array CNV types for {tool2_label}\n")
                    for acnvtype in comparisondata[comparepoint]["tool2"]:
                        outfile.write(acnvtype + ": " + str(comparisondata[comparepoint]["tool2"][acnvtype]) + "\n")
                    # outfile.write("\n\n")

                # Write data for the overlap in shared array CNVs
                # if comparepoint == "Shared overlap":
        file_written = True
    except IOError:
        print("Could not write comparison data to output file")
    finally:
        return file_written


def write_fp_comparison(outfilepath, comparisondata, tool1_label, tool2_label):
    try:
        file_written = False
        comparison_points = ["Total fps", "Shared fps", "Overlapping fps", "Unique fps"]

        with open(outfilepath, 'w') as outfile:
            outfile.write(f"Results for comparison between {tool1_label} and {tool2_label}\n\n")

            for comparepoint in comparison_points:
                if comparepoint == "Total fps":
                    outfile.write("[-Total False Positive calls-]\n")
                    outfile.write(f"{tool1_label}: {comparisondata[comparepoint][0]}\n")
                    outfile.write(f"{tool2_label}: {comparisondata[comparepoint][1]}\n")
                    outfile.write("\n\n")

                if comparepoint == "Shared fps":
                    num_of_shared_fps = 0
                    outfile.write("[-Shared False Positive calls-]\n")
                    for samplename in comparisondata[comparepoint]:
                        shared_sample_fps = len(comparisondata[comparepoint][samplename])
                        num_of_shared_fps += shared_sample_fps
                        outfile.write("* [" + samplename + "]: {" +str(shared_sample_fps)+ "}\t" + ", ".join(comparisondata[comparepoint][samplename]) + "\n")
                    outfile.write(f"{tool1_label} and {tool2_label} both identified {num_of_shared_fps} False Positive calls\n")
                    outfile.write("\n\n")

                if comparepoint == "Overlapping fps":
                    num_of_overlaps = 0
                    outfile.write("[-Overlapping False Positive calls-]\n")
                    for samplename in comparisondata[comparepoint]:
                        sample_overlaps = len(comparisondata[comparepoint][samplename])
                        num_of_overlaps += sample_overlaps
                        outfile.write("* [" +samplename+ "]: {" +str(sample_overlaps)+ "}\t")
                        for overlappingfps in comparisondata[comparepoint][samplename]:
                            outfile.write("(" +overlappingfps[0]+ " - " +overlappingfps[1]+ "), ")
                        outfile.write("\n")
                    outfile.write(f"{tool1_label} and {tool2_label} have {num_of_overlaps} overlapping False Positive calls\n")
                    outfile.write("\n\n")

                if comparepoint == "Unique fps":
                    outfile.write(f"[-Unique False Positive calls found by {tool1_label}-]\n")
                    num_of_tool1_fps = 0
                    for samplename in comparisondata[comparepoint]["tool1"]:
                        tool1_sample_ufps = len(comparisondata[comparepoint]["tool1"][samplename])
                        num_of_tool1_fps += tool1_sample_ufps
                        outfile.write("* [" +samplename+ "]: {" +str(tool1_sample_ufps)+ "}\t" + ", ".join(comparisondata[comparepoint]["tool1"][samplename]) + "\n")
                    outfile.write(f"{tool1_label} found {num_of_tool1_fps} unique False Positive calls\n\n")

                    outfile.write(f"[-Unique False Positive calls found by {tool2_label}-]\n")
                    num_of_tool2_fps = 0
                    for samplename in comparisondata[comparepoint]["tool2"]:
                        tool2_sample_ufps = len(comparisondata[comparepoint]["tool2"][samplename])
                        num_of_tool2_fps += tool2_sample_ufps
                        outfile.write("* [" +samplename+ "]: {" +str(tool2_sample_ufps)+ "}\t" + ", ".join(comparisondata[comparepoint]["tool2"][samplename]) + "\n")
                    outfile.write(f"{tool2_label} found {num_of_tool2_fps} unique False Positive calls\n\n")
        file_written = True
    except IOError:
        print("Could not write false positive comparison data to output file")
    finally:
        return file_written


def write_tp_comparison(outfilepath, comparisondata, tool1_label, tool2_label):
    try:
        file_written = False
        comparison_points = ["Total tps", "Shared tps", "Overlapping tps", "Unique tps"]

        with open(outfilepath, 'w') as outfile:
            outfile.write(f"Results for comparison between {tool1_label} and {tool2_label}\n\n")

            for comparepoint in comparison_points:
                if comparepoint == "Total tps":
                    outfile.write("[-Total True Positive calls-]\n")
                    outfile.write(f"{tool1_label}: {comparisondata[comparepoint][0]}\n")
                    outfile.write(f"{tool2_label}: {comparisondata[comparepoint][1]}\n")
                    outfile.write("\n\n")

                if comparepoint == "Shared tps":
                    num_of_shared_fps = 0
                    outfile.write("[-Shared True Positive calls-]\n")
                    for samplename in comparisondata[comparepoint]:
                        shared_sample_fps = len(comparisondata[comparepoint][samplename])
                        num_of_shared_fps += shared_sample_fps
                        outfile.write("* [" + samplename + "]: {" +str(shared_sample_fps)+ "}\t" + ", ".join(comparisondata[comparepoint][samplename]) + "\n")
                    outfile.write(f"{tool1_label} and {tool2_label} both identified {num_of_shared_fps} True Positive calls\n")
                    outfile.write("\n\n")

                if comparepoint == "Overlapping tps":
                    num_of_overlaps = 0
                    outfile.write("[-Overlapping True Positive calls-]\n")
                    for samplename in comparisondata[comparepoint]:
                        sample_overlaps = len(comparisondata[comparepoint][samplename])
                        num_of_overlaps += sample_overlaps
                        outfile.write("* [" +samplename+ "]: {" +str(sample_overlaps)+ "}\t")
                        for overlappingfps in comparisondata[comparepoint][samplename]:
                            outfile.write("(" +overlappingfps[0]+ " - " +overlappingfps[1]+ "), ")
                        outfile.write("\n")
                    outfile.write(f"{tool1_label} and {tool2_label} have {num_of_overlaps} overlapping True Positive calls\n")
                    outfile.write("\n\n")

                if comparepoint == "Unique tps":
                    outfile.write(f"[-Unique True Positive calls found by {tool1_label}-]\n")
                    num_of_tool1_fps = 0
                    for samplename in comparisondata[comparepoint]["tool1"]:
                        tool1_sample_ufps = len(comparisondata[comparepoint]["tool1"][samplename])
                        num_of_tool1_fps += tool1_sample_ufps
                        outfile.write("* [" +samplename+ "]: {" +str(tool1_sample_ufps)+ "}\t" + ", ".join(comparisondata[comparepoint]["tool1"][samplename]) + "\n")
                    outfile.write(f"{tool1_label} found {num_of_tool1_fps} unique True Positive calls\n\n")

                    outfile.write(f"[-Unique True Positive calls found by {tool2_label}-]\n")
                    num_of_tool2_fps = 0
                    for samplename in comparisondata[comparepoint]["tool2"]:
                        tool2_sample_ufps = len(comparisondata[comparepoint]["tool2"][samplename])
                        num_of_tool2_fps += tool2_sample_ufps
                        outfile.write("* [" +samplename+ "]: {" +str(tool2_sample_ufps)+ "}\t" + ", ".join(comparisondata[comparepoint]["tool2"][samplename]) + "\n")
                    outfile.write(f"{tool2_label} found {num_of_tool2_fps} unique True Positive calls\n\n")
        file_written = True
    except IOError:
        print("Could not write false positive comparison data to output file")
    finally:
        return file_written
