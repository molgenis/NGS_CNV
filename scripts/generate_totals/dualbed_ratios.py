from classes.gatkcall import GatkCall

def read_dualbed_data(infileloc):
    """Read and return dualBED classified data.

    Parameters
    ----------
    dualbedfile : str
        Path to dualBED classified file

    Returns
    -------
    dualbed_data : dict
        dualBED classified data saved per sample, per dualBED label
    """
    dualbed_data = {}
    dualbed_data["Shared"] = {}
    dualbed_data["Overlapping"] = {}
    dualbed_data["Unique"] = {}

    try:
        with open(infileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")

                gatkreg = get_cnv_region(filelinedata[1])
                gatkgenes = filelinedata[13].split(":")
                gatkcall = GatkCall(filelinedata[0], gatkreg[0], gatkreg[1], gatkreg[2], filelinedata[2], int(filelinedata[3]), filelinedata[4], filelinedata[9], filelinedata[10], int(filelinedata[11]), int(filelinedata[12]), gatkgenes, fileline.strip())

                if filelinedata[-1] == "Shared":
                    if filelinedata[0] not in dualbed_data["Shared"]:
                        dualbed_data["Shared"][filelinedata[0]] = []
                    dualbed_data["Shared"][filelinedata[0]].append(gatkcall)

                elif filelinedata[-1] == "Overlapping":
                    if filelinedata[0] not in dualbed_data["Overlapping"]:
                        dualbed_data["Overlapping"][filelinedata[0]] = []
                    dualbed_data["Overlapping"][filelinedata[0]].append(gatkcall)

                else:
                    if filelinedata[0] not in dualbed_data["Unique"]:
                        dualbed_data["Unique"][filelinedata[0]] = []
                    dualbed_data["Unique"][filelinedata[0]].append(gatkcall)
    except IOError:
        print("Could not read input file")
    finally:
        return dualbed_data


def merge_overlapping(normaloverlaps, hcoverlaps):
    """Merge normal and HighConfident overlapping calls into a single set.

    Parameters
    ----------
    normaloverlaps : dict
        Normal dualBED overlapping calls
    hcoverlaps : dict
        HighConfident dualBED overlapping calls

    Returns
    -------
    overlapping_calls : dict
        Merged normal and HighConfident dualBED overlapping calls
    """
    overlapping_calls = {}

    # Iterate over the normal overlapping dualBED calls
    for samplename in normaloverlaps:
        if samplename not in overlapping_calls:
            overlapping_calls[samplename] = []
        for nocnv in normaloverlaps[samplename]:
            overlapping_calls[samplename].append(nocnv)

    # Iterate over the high confident overlapping dualBED calls
    for samplename in hcoverlaps:
        if samplename not in overlapping_calls:
            overlapping_calls[samplename] = []
        for hcocnv in hcoverlaps[samplename]:
            overlapping_calls[samplename].append(hcocnv)
    return overlapping_calls


def form_shared_filter(shareddata):
    """Create and return shared array CNV filter for generating ratio totals for overlapping calls

    Parameters
    ----------
    shareddata : dict
        CNV calls with the dualBED label Shared

    Returns
    -------
    shared_filter
        List of all shared array CNV
    """
    shared_filter = {}
    for samplename in shareddata:
        for scnv in shareddata[samplename]:
            if scnv.classification == "POSITIVE" or scnv.classification == "TRUE POSITIVE":
                if samplename not in shared_filter:
                    shared_filter[samplename] = []
                shared_filter[samplename].append(scnv.arraycnv)
    return shared_filter


def generate_overlapping_totals(overlappingcalls, sharedfilter, tpperacnv):
    """Generate and return classification totals for overlapping calls

    Parameters
    ----------
    overlappingcalls
    sharedfilter
    tpperacnv : bool
        Count TPs only per array CNV if true

    Returns
    -------
    overlapping_totals : dict
        Overlapping totals
    """
    overlapping_totals = {}
    arraycnvs_found = {}
    for samplename in overlappingcalls:
        for ocnv in overlappingcalls[samplename]:
            classlabel = determine_clasification_totals_label(ocnv.classification)
            count_ocnv = True
            if classlabel == "True Positive":
                # Check if this TP was already found in Shared
                if samplename in sharedfilter:
                    if ocnv.arraycnv in sharedfilter[samplename]:
                        count_ocnv = False

            # Check if the overlapping CNV should be counted
            if count_ocnv:
                if classlabel not in overlapping_totals:
                    overlapping_totals[classlabel] = 0

                if classlabel == "True Positive" and tpperacnv:
                    if samplename not in arraycnvs_found:
                        arraycnvs_found[samplename] = []
                    if ocnv.arraycnv not in arraycnvs_found[samplename]:
                        overlapping_totals[classlabel] += 1
                        arraycnvs_found[samplename].append(ocnv.arraycnv)
                else:
                    overlapping_totals[classlabel] += 1
    return overlapping_totals


def determine_clasification_totals_label(classification):
    """Determine the classification label to assign

    Parameters
    ----------
    classification : str
        Assigned classification by classification.py

    Returns
    -------
    classlabel : str
        Label to use for the totals
    """
    classlabel = ""
    if "ARRAY NON-" in classification:
        classlabel = "Array Non-Informative"
    if "WES NON-" in classification:
        classlabel = "WES Non-Informative"
    if "ARRAY & WES NON-" in classification:
        classlabel = "Array & WES Non-Informative"
    if classification == "FALSE POSITIVE":
        classlabel = "False Positive"
    if classification == "POSITIVE":
        classlabel = "True Positive"
    if classification == "TRUE POSITIVE":
        classlabel = "True Positive"
    return classlabel


def write_dualbed_ratios(outfileloc, sharedratios, overlappingratios, nuniqueratios, hcuniqueratios):
    """Write dualBED ratios to an output file.

    Parameters
    ----------
    outfileloc : str
        Path to write output file to
    sharedratios : dict
        Ratios for the Shared dualBED CNV calls
    overlappingratios : dict
        Ratios for the Overlapping dualBED CNV calls
    nuniqueratios : dict
        Ratios for the Unique normal dualBED CNV calls
    hcuniqueratios : dict
        Ratios for the Unique HighConfident dualBED CNV calls

    Returns
    -------
    file_written : bool
        True if output file has succesfully been written, False if not
    """
    file_written = False
    ratiolabels = ["True Positive", "False Positive", "Array Non-Informative", "WES Non-Informative", "Array & WES Non-Informative"]
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("[-Shared-]\n")
            for ratiolabel in ratiolabels:
                if ratiolabel in sharedratios:
                    outfile.write(f"{ratiolabel}: {sharedratios[ratiolabel]}\n")
            outfile.write("\n")

            # Write the Overlapping ratios to file
            outfile.write("[-Overlapping-]\n")
            for ratiolabel in ratiolabels:
                if ratiolabel in overlappingratios:
                    outfile.write(f"{ratiolabel}: {overlappingratios[ratiolabel]}\n")
            outfile.write("\n")

            # Write the normal Unique ratios to file
            outfile.write("[-Unique Normal-]\n")
            for ratiolabel in ratiolabels:
                if ratiolabel in nuniqueratios:
                    outfile.write(f"{ratiolabel}: {nuniqueratios[ratiolabel]}\n")
            outfile.write("\n")

            # Write the high confident Unique ratios to file
            outfile.write("[-Unique HighConfident-]\n")
            for ratiolabel in ratiolabels:
                if ratiolabel in hcuniqueratios:
                    outfile.write(f"{ratiolabel}: {hcuniqueratios[ratiolabel]}\n")
        file_written = True
    except IOError:
        print("Could not write ratio data to file")
    finally:
        return file_written


def get_cnv_region(arraycnvregion):
    """Split and return the array CNV region into the chrom, start and stop.

    Parameters
    ----------
    arraycnvregion : str
        Array CNV region data

    Returns
    -------
    list of str and int
        Array CNV region data
    """
    arrayregion = arraycnvregion.replace(",", "")
    arrayregiondata = arrayregion.split(":")
    arraystart = arrayregiondata[1].split("-")[0]
    arrayend = arrayregiondata[1].split("-")[1]
    return [arrayregiondata[0], arraystart, arrayend]
