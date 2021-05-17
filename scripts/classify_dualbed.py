import argparse
import utils.filereaders as ufr


def get_params():
    dualbed_args = argparse.ArgumentParser()
    dualbed_args.add_argument("-1", "--infile1", dest="infile1", help="Path to classified CNV calls for the normal BED file")
    dualbed_args.add_argument("-2", "--infile2", dest="infile2", help="Path to classified CNV calls for the High Confident BED file")
    dualbed_args.add_argument("-o", "--outdir", dest="outdir", help="Path to output directory to write output files to")
    return vars(dualbed_args.parse_args())


def get_call_regions(classificationdata):
    """Return all calls as regions strings per sample

    Parameters
    ----------
    classificationdata  : dict
        Classification data of CNV calling analysis

    Returns
    -------
    classification_regions : dict
        Classified calls as strings per sample
    """
    classification_regions = {}
    for samplename in classificationdata:
        for cnvcall in classificationdata[samplename]:
            if samplename not in classification_regions:
                classification_regions[samplename] = []
            classification_regions[samplename].append(cnvcall.get_region_str())
    return classification_regions


def get_shared_calls(normaldata, highconfidentdata):
    """Determine calls that are exactly the same for normal and HC bed file.

    Parameters
    ----------
    normaldata : dict
    highconfidentdata : dict
    """
    shared_calls = {}
    for samplename in set(normaldata.keys()) & set(highconfidentdata.keys()):
        shared_sample_calls = set(normaldata[samplename]) & set(highconfidentdata[samplename])
        if len(shared_sample_calls) > 0:
            shared_calls[samplename] = shared_sample_calls
    return shared_calls


def get_overlapping_calls(normaldata, highconfidentdata, sharedcalls):
    """Deterine calls in the normal and HC bed file that overlap.

    Parameters
    ----------
    normaldata : dict
    highconfidentdata : dict
    sharedcalls : dict
        Shared calls between normal and HC

    Returns
    -------
    overlapping_calls : dict
        Overlapping calls between normal and HC per sample
    """
    overlapping_calls = {}
    for samplename in set(normaldata.keys()) & set(highconfidentdata.keys()):
        for cnvcall in normaldata[samplename]:
            if samplename not in sharedcalls:
                overlapregions = regions_overlap(cnvcall.get_region_str(), highconfidentdata[samplename])
                if len(overlapregions) > 0:
                    if samplename not in overlapping_calls:
                        overlapping_calls[samplename] = []
                    overlapping_calls[samplename].extend(overlapregions)
            else:
                if cnvcall.get_region_str() not in sharedcalls[samplename]:
                    overlapregions = regions_overlap(cnvcall.get_region_str(), highconfidentdata[samplename])
                    if len(overlapregions) > 0:
                        if samplename not in overlapping_calls:
                            overlapping_calls[samplename] = []
                        overlapping_calls[samplename].extend(overlapregions)
    return overlapping_calls


def regions_overlap(normalregion, hccalls):
    overlapping_calls = []
    normal_c = normalregion.split(":")[0]
    normal_s = int(normalregion.split(":")[1].split("-")[0])
    normal_e = int(normalregion.split(":")[1].split("-")[1])

    for hccall in hccalls:
        hcregion = hccall.get_region_str()
        hc_c = hcregion.split(":")[0]
        if normal_c == hc_c:
            hc_s = int(hcregion.split(":")[1].split("-")[0])
            hc_e = int(hcregion.split(":")[1].split("-")[1])
            if normal_s <= hc_e and hc_s <= normal_e:
                overlapping_calls.append((normalregion, hcregion))
    return overlapping_calls


def get_unique_calls(normaldata, highconfidentdata, overlappingcalls):
    unique_calls = {}
    unique_calls["normal"] = {}
    unique_calls["hc"] = {}
    
    for samplename in set(normaldata.keys()) | set(highconfidentdata.keys()):
        overlapfilter = []
        if samplename in overlappingcalls:
            overlapfilter = [x[0] for x in overlappingcalls[samplename]]
            hcoverlaps = [x[1] for x in overlappingcalls[samplename]]
            overlapfilter.extend(hcoverlaps)

        normal_regions = [cnvcall.get_region_str() for cnvcall in normaldata[samplename]]
        hc_regions = [cnvcall.get_region_str() for cnvcall in highconfidentcall[samplename]]

        normal_unique_calls = set(normal_regions) - set(hc_regions)
        hc_unique_calls = set(hc_regions) - set(normal_regions)
        normal_unique_calls = normal_uniqe_calls - set(overlapfilter)
        hc_unique_calls = hc_unique_calls - set(overlapfilter)

        if len(normal_unique_calls) > 0:
            unique_calls["normal"][samplename] = normal_unique_calls
        if len(hc_unique_calls) > 0:
            unique_calls["hc"][samplename] = hc_unique_calls
    return unique_calls


def write_dualbed_file(outfilepath, classificationdata, sharedcalls, overlappingcalls, uniquecalls):
    """Write classification output file with added dualBED label for each call.

    Parameters
    ----------
    outfilepath : str
        Path to write output file to
    classificationdata : dict
        Classification data to write to file with added label
    sharedcalls : dict
        Calls shared with the other BED file
    overlappingcalls : dict
        Calls overlapping with the other BED file
    uniquecalls : dict
        Calls unique to this BED file
    """
    file_written = False
    try:
        with open(outfilepath, 'w') as outfile:
            outfile.write("Sample\tGATK4_CNV\tGATK4_Call\tGATK4_Size\tArray_CNV\tArray_Call\tArray_Size\tHangover_L\tHangover_R\tCall_Result\tClassification\t#_Exons\t#_Probes\tGATK4_genes\tArray_genes\tGATK4_UGenes\tArray_UGenes\tDualBED\n")
            for samplename in classificationdata:
                shared_sample_calls = set()
                overlapping_sample_calls = set()
                unique_sample_calls = set()

                if samplename in sharedcalls:
                    shared_sample_calls = set(sharedcalls[samplename])
                if samplename in overlappingcalls:
                    overlapping_sample_calls = [x[0] for x in overlappingcalls[samplename]]
                    hc_sample_overlap = [x[1] for x in overlappingcalls[samplename]]
                    overlapping_sample_calls.extend(hc_sample_overlap)
                if samplename in uniquecalls:
                    unique_sample_calls = set(uniquecalls[samplename])

                for cnvcall in classificationdata[samplename]:
                    dualbedlabel = get_dualbed_label(cnvcall.get_region_str(), shared_sample_calls, overlapping_sample_calls, unique_sample_calls)
                    outfile.write(f"{cnvcall.call_line}\t{dualbedlabel}\n")
        file_written = True
    except IOError:
        print("Could not write output file")
    finally:
        return file_written


def get_dualbed_label(cnvcallregion, sharedcalls, overlappingcalls, uniquecalls):
    """Determine and return the dualbed label for a specified CNV call.

    Parameters
    ----------
    cnvcallregion : str
        CNV call region to get dualbed label for
    sharedcalls : list of str
        List of shared CNV calls as region strings for a single sample
    overlappingcalls : lst of str
        List of overlapping CNV calls as region strings for a single sample
    uniquecalls : list of str
        List of unique CNV calls as region strings for a single sample
    """
    if cnvcallregion in sharedcalls:
        return "Shared"
    elif cnvcallregion in overlappingcalls:
        return "Overlapping"
    elif cnvcallregion in uniquecalls:
        return "Unique"


def main():
    """Perform the work"""
    dualbed_params = get_params()
    normal_data = ufr.read_classification_file(dualbed_params["infile1"])
    highconfident_data = ufr.read_classification_file(dualbed_params["infile2"])

    # Obtain the normal and high confident calls as region strings
    normal_regions = get_call_regions(normal_data)
    hc_regions = get_calls_regions(highconfident_data)

    # Determine the shared, overlapping and unique calls between the normal and high confident classification data
    shared_calls = get_shared_calls(normal_regions, hc_regions)
    overlapping_calls = get_overlapping_calls(normal_regions, hc_regions, shared_calls)
    unique_calls = get_unique_calls(normal_regions, hc_regions, overlapping_calls)

    # Write the normal classification data with the added dualbed label (Shared, Overlapping, Unique)
    normal_outpath = dualbed_params["outdir"]+ "/" +dualbed_params["infile1"].split("/")[-1]
    write_dualbed_file(normal_outpath, normal_regions, shared_calls, overlapping_calls, unique_calls["normal"])

    # Write the High Confident classification data with the added dualbed label (Shared, Overlapping, Unique)
    highconfident_outpath = dualbed_params["outdir"]+ "/" +dualbed_params["infile2"].split("/")[-1]
    write_dualbed_file(highconfident_outpath, hc_regions, shared_calls, overlapping_calls, unique_calls["hc"])


if __namne__ == "__main__":
    main()
