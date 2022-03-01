#!/usr/bin/env python
def determine_dualbed_acnvs_found(infileloc, arraycnvs, dualbedlabel):
    """Determine and return which array CNVs have been found.

    Parameters
    ----------
    infileloc : str
        Path to classifications file
    arraycnvs : dict
        Read array CNV data from array CNV file

    Returns
    -------
    array_cnv_found : dict
        Found array CNV regions saved per sample
    """
    already_found = {}
    array_cnvs_found = {}
    try:
        with open(infileloc, 'r') as infile:
            next(infile)
            for inline in infile:
                inlinedata = inline.strip().split("\t")

                if inlinedata[4] != "NA":
                    if inlinedata[-1] == dualbedlabel:
                        if inlinedata[0] not in array_cnvs_found:
                            array_cnvs_found[inlinedata[0]] = []

                        if inlinedata[0] not in already_found:
                            already_found[inlinedata[0]] = []

                        if inlinedata[4] not in already_found[inlinedata[0]]:
                            array_cnvs_found[inlinedata[0]].append(arraycnvs[inlinedata[0]][inlinedata[4]])

                        if inlinedata[4] not in already_found[inlinedata[0]]:
                            already_found[inlinedata[0]].append(inlinedata[4])
    except IOError:
        print(f"Could not open {infileloc}")
    finally:
        return array_cnvs_found


def determine_dualbed_acnvs_found_2(infileloc, arraycnvs, dualbedlabel, soufilter):
    already_found = {}
    array_cnvs_found = {}
    try:
        with open(infileloc, 'r') as infile:
            next(infile)
            for inline in infile:
                inlinedata = inline.strip().split("\t")

                if inlinedata[4] != "NA":
                    if inlinedata[-1] == dualbedlabel:
                        # Check if the current array CNV is in the SOU filter
                        if not acnv_in_soufilter(inlinedata[0], inlinedata[4], soufilter):
                            if inlinedata[0] not in array_cnvs_found:
                                array_cnvs_found[inlinedata[0]] = []

                            if inlinedata[0] not in already_found:
                                already_found[inlinedata[0]] = []

                            if inlinedata[4] not in already_found[inlinedata[0]]:
                                array_cnvs_found[inlinedata[0]].append(arraycnvs[inlinedata[0]][inlinedata[4]])

                            if inlinedata[4] not in already_found[inlinedata[0]]:
                                already_found[inlinedata[0]].append(inlinedata[4])
    except IOError:
        print(f"Could not open {infileloc}")
    finally:
        return array_cnvs_found


def acnv_in_soufilter(samplename, acnvregion, soufilter):
    """Determine whether an array CNV is in the SOU filter.

    Parameters
    ----------
    samplename : str
        Name of the sample
    acnvregion : str
        Array CNV region as string representation
    soufilter : dict
        Shared, Overlapping, Unique filter

    Returns
    -------
    acnv_in_filter : bool
        True if array CNV is in SOU filter, False if not
    """
    acnv_in_filter = False
    if samplename in soufilter:
        acnv_in_filter = acnvregion in soufilter[samplename]
    return acnv_in_filter


def summarize_found_arraycnv_types(sharedacnvs, overlappingacnvs, uniqueacnvs):
    typecounts = {}
    typecounts["Shared"] = {}
    typecounts["Overlapping"] = {}
    typecounts["Unique"] = {}
    acnvs_processed = {}

    # Summarize Shared
    for samplename in sharedacnvs:
        for sfacnv in sharedacnvs[samplename]:
            if sfacnv.cnv_class not in typecounts["Shared"]:
                typecounts["Shared"][sfacnv.cnv_class] = 0

            if samplename not in acnvs_processed:
                acnvs_processed[samplename] = []

            if sfacnv.get_region() not in acnvs_processed[samplename]:
                typecounts["Shared"][sfacnv.cnv_class] = typecounts["Shared"][sfacnv.cnv_class] + 1
                acnvs_processed[samplename].append(sfacnv.get_region())

    # Summarize Overlapping
    for samplename in overlappingacnvs:
        for ofacnv in overlappingacnvs[samplename]:
            if ofacnv.cnv_class not in typecounts["Overlapping"]:
                typecounts["Overlapping"][ofacnv.cnv_class] = 0

            if samplename not in acnvs_processed:
                acnvs_processed[samplename] = []

            if ofacnv.get_region() not in acnvs_processed:
                typecounts["Overlapping"][ofacnv.cnv_class] = typecounts["Overlapping"][ofacnv.cnv_class] + 1
                acnvs_processed[samplename].append(ofacnv.get_region())

    # Summarize Unique
    for samplename in uniqueacnvs:
        for ufacnv in uniqueacnvs[samplename]:
            if ufacnv.cnv_class not in typecounts["Unique"]:
                typecounts["Unique"][ufacnv.cnv_class] = 0

            if samplename not in acnvs_processed:
                acnvs_processed[samplename] = []

            if ufacnv.get_region() not in acnvs_processed[samplename]:
                typecounts["Unique"][ufacnv.cnv_class] = typecounts["Unique"][ufacnv.cnv_class] + 1
                acnvs_processed[samplename].append(ufacnv.get_region())
    return typecounts


def determine_dualbed_acnvs_missed(arraycnvs, sfound_data, ofound_data, ufound_data):
    """Make and return a list of all the missed array CNVs.

    Parameters
    ----------
    arraycnvs : dict
        Array CNVs
    sfound_data : dict
        All shared array CNVs found
    ofound_data : dict
        All overlapping array CNVs found
    ufound_data : dict
        All unique array CNVs found

    Returns
    -------
    missed_array_cnvs : dict
        Missed array CNVs per sample
    """
    missed_array_cnvs = {}
    for samplename in arraycnvs:
        for acnv in arraycnvs[samplename]:
            sfound = acnv_in_soufilter(samplename, arraycnvs[samplename][acnv].get_region(), sfound_data)
            ofound = acnv_in_soufilter(samplename, arraycnvs[samplename][acnv].get_region(), ofound_data)
            ufound = acnv_in_soufilter(samplename, arraycnvs[samplename][acnv].get_region(), ufound_data)

            if not sfound and not ofound and not ufound:
                if samplename not in missed_array_cnvs:
                    missed_array_cnvs[samplename] = []
                missed_array_cnvs[samplename].append(arraycnvs[samplename][acnv])
    return missed_array_cnvs


def summarize_missed_arraycnv_types(missed_acnvs):
    typecounts = {}
    for samplename in missed_acnvs:
        for acnv in missed_acnvs[samplename]:
            if acnv.cnv_class not in typecounts:
                typecounts[acnv.cnv_class] = 0
            typecounts[acnv.cnv_class] = typecounts[acnv.cnv_class] + 1
    return typecounts


def write_foundmissed_acnvs(missedfound_arraycnvs, outfileloc):
    """Write found/missed array cnvs to file

    Parameters
    ----------
    missedfound_arraycnvs : dict
        Array CNVs to write
    outfileloc : str
        Path to write output file to
    """
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


def write_found_summary(found_summary, outfileloc):
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("[-Shared-]\n")
            outfile.write("Cnv_Type\tCount\n")
            for scnvtype in found_summary["Shared"]:
                typecount = found_summary["Shared"][scnvtype]
                outfile.write(f"{scnvtype}\t{typecount}\n")
            outfile.write("\n")

            outfile.write("[-Overlapping-]\n")
            outfile.write("Cnv_Type\tCount\n")
            for ocnvtype in found_summary["Overlapping"]:
                typecount = found_summary["Overlapping"][ocnvtype]
                outfile.write(f"{ocnvtype}\t{typecount}\n")
            outfile.write("\n")

            outfile.write("[-Unique-]\n")
            outfile.write("Cnv_Type\tCount\n")
            for ucnvtype in found_summary["Unique"]:
                typecount = found_summary["Unique"][ucnvtype]
                outfile.write(f"{ucnvtype}\t{typecount}\n")
        file_written = True
    except IOError:
        print("")
    finally:
        return file_written


def write_missed_summary(missed_summary, outfileloc):
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Cnv_Type\tCount\n")
            for mcnvtype in missed_summary:
                outfile.write(f"{mcnvtype}\t{missed_summary[mcnvtype]}\n")
    except IOError:
        print("")
    finally:
        return file_written


def make_sou_found_filter(sou_data):
    """Make and return a list of already found Shared, Overlapping or Unique array CNVs to use as a filter.

    Parameters
    ----------
    sou_data : dict
        Shared, Overlapping or Unique found array CNVs

    Returns
    -------
    sou_filter : dict
        Filter to use
    """
    sou_filter = {}
    for samplename in sou_data:
        if samplename not in sou_filter:
            sou_filter[samplename] = []
        for acnv in sou_data[samplename]:
            sou_filter[samplename].append(acnv.get_region())
    return sou_filter

