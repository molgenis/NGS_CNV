# Import generate totals scripts
import generate_totals.classifications as gtc
import generate_totals.array_cnvs as gtac


def perform_comparison(tool1_label, tool1_data, tool2_label, tool2_data, arraydata):
    """Compare classified CNV calls from two tools.

    Parameters
    ----------
    tool1_label : str
        Name for the first tool
    tool1_data : dict
        Classified CNV calls for the first tool
    tool2_label : str
        Name for the second tool
    tool2_data : dict
        Classified CNV calls for the second tool

    Returns
    -------
    comparison_data : dict
        Comparison summary data
    """
    comparison_data = {}

    print("COMPARISON: [Compare total calls made]")
    total_calls = compare_total_calls(tool1_data, tool2_data)
    comparison_data["Total calls"] = total_calls

    print("COMPARISON: [Compare total calls per classifcation label]")
    call_types = compare_call_types(tool1_data, tool2_data)
    comparison_data["Call types"] = call_types

    print("COMPARISON: [Determine the number of found array CNVs]")
    found_arraycnvs = compare_found_array_cnvs(tool1_data, tool2_data)
    comparison_data["Found array CNVs"] = found_arraycnvs

    print("COMPARISON: [Determine the number of shared arryay CNVs]")
    shared_arraycnvs = determine_shared_array_cnvs(found_arraycnvs)
    comparison_data["Shared array CNVs"] = shared_arraycnvs

    print("[Determine the shared array CNV types]")
    shared_acnv_types = determine_shared_acnv_types(shared_arraycnvs, arraydata)
    comparison_data["Shared array types"] = shared_acnv_types

    print("COMPARISON: [Determine array CNVs found by only one tool but not the other]")
    unique_arraycnvs = determine_unique_array_cnvs(found_arraycnvs)
    comparison_data["Unique array CNVs"] = unique_arraycnvs

    print("[Determine the type of uniquely found array CNV types per tool]")
    unique_acnv_types = determine_unique_acnv_types(unique_arraycnvs, arraydata)
    comparison_data["Unique array types"] = unique_acnv_types

    # print("COMPARISON: [Determine the amount of overlap with shared array CNVs]")
    return comparison_data


def compare_total_calls(tool1_data, tool2_data):
    """Compare and return the total number of calls for two tools.

    Parameters
    ----------
    tool1_data : dict
        Classified CNV calls for the first tool
    tool2_data : dict
        Classified CNV calls for the second tool

    Returns
    -------
    list of int
        Total calls for tool 1 and 2
    """
    tool1_callnum = determine_total_calls(tool1_data)
    tool2_callnum = determine_total_calls(tool2_data)
    return [tool1_callnum, tool2_callnum]
    # return [determine_total_calls(tool1_data), determine_total_calls(tool2_data)]


def determine_total_calls(tooldata):
    """Determine and return the total number of calls.

    Parameters
    ----------
    tooldata : dict
        Classified CNV calls

    Returns
    -------
    tool_callnum : int
        Total number of CNV calls made.
    """
    tool_callnum = 0
    for samplename in tooldata:
        tool_callnum += len(tooldata[samplename])
    return tool_callnum


def compare_call_types(tool1_data, tool2_data):
    """Compares the number of TP, FP, ANI, WNI, AWNI calls between the two tools.

    Parameters
    ----------
    tool1_data : dict
        Classified CNV calls for the first tool
    tool2_data : dict
        Classified CNV calls for the second tool

    Returns
    -------
    classification_numbers : dict
        Number of calls for TP, FP, ANI, WNI, AWNI
    """
    classification_labels = ["True Positive", "False Positive", "Array Non-Informative", "WES Non-Informative", "Array & WES Non-Informative"]
    tool1_classification_numbers = gtc.generate_classification_totals(tool1_data)
    tool2_classification_numbers = gtc.generate_classification_totals(tool2_data)

    classification_numbers = {}
    for classlabel in classification_labels:
        classification_numbers[classlabel] = [0, 0]
        if classlabel in tool1_classification_numbers:
            classification_numbers[classlabel][0] = tool1_classification_numbers[classlabel]
        if classlabel in tool2_classification_numbers:
            classification_numbers[classlabel][1] = tool2_classification_numbers[classlabel]
    return classification_numbers


def compare_found_array_cnvs(tool1_data, tool2_data):
    """Compare the number of found array CNVs.

    Parameters
    ----------
    tool1_data : dict
        Classified CNV calls for the first tool
    tool2_data : dict
        Classified CNV calls for the second tool

    Returns
    -------
    found_array_cnvs : dict
        Array CNVs found for both tools
    """
    tool1_arraycnvs = gtac.array_cnvs_found(tool1_data)
    tool2_arraycnvs = gtac.array_cnvs_found(tool2_data)
    found_array_cnvs = {"tool1": tool1_arraycnvs,
                        "tool2": tool2_arraycnvs}
    return found_array_cnvs


def determine_shared_array_cnvs(found_arraycnv_data):
    """Determine which array CNVs are shared between the two tools.

    Parameters
    ----------
    found_arraycnv_data : dict
        Found array CNVs for both tools

    Returns
    -------
    shared_array_cnvs : dict
        Found array CNVs shared between the two tools
    """
    shared_array_cnvs = {}
    for samplename in set(found_arraycnv_data["tool1"].keys()) & set(found_arraycnv_data["tool2"].keys()):
        shared_array_cnvs[samplename] = set(found_arraycnv_data["tool1"][samplename]) & set(found_arraycnv_data["tool2"][samplename])
    return shared_array_cnvs


def determine_shared_acnv_types(shared_arraycnvs, arraycnvdata):
    """Determine the types of the shared array CNVs.

    Parameters
    ----------
    shared_arraycnvs : dict
        Array CNVs found by both tools
    arraycnvdata : dict
        Array CNV data

    Returns
    -------
    shared_acnv_types : dict
        Count per array CNV type
    """
    shared_acnv_types = {}
    for samplename in shared_arraycnvs:
        for shared_acnv in shared_arraycnvs[samplename]:
            if shared_acnv in arraycnvdata[samplename]:
                array_type = arraycnvdata[samplename][shared_acnv].cnv_class
                if array_type not in shared_acnv_types:
                    shared_acnv_types[array_type] = 0
                shared_acnv_types[array_type] += 1
    return shared_acnv_types


def determine_unique_array_cnvs(found_array_cnv_data):
    """Determine which array CNVs are unique for each tool.

    Parameters
    ----------
    found_array_cnv_data : dict
        Array CNVs found by both tools

    Returns
    -------
    unique_array_cnvs : dict
        Found array CNVs unique to each tool
    """
    unique_array_cnvs = {"tool1": {}, "tool2": {}}
    for samplename in set(found_array_cnv_data["tool1"].keys()) | set(found_array_cnv_data["tool2"].keys()):
        tool1_acnvs, tool2_acnvs = set(), set()

        # Gather the array CNVs found for the current sample.
        if samplename in found_array_cnv_data["tool1"]:
            tool1_acnvs = set(found_array_cnv_data["tool1"][samplename])
        if samplename in found_array_cnv_data["tool2"]:
            tool2_acnvs = set(found_array_cnv_data["tool2"][samplename])

        # Determine the unique array CNVs for both tools
        tool1_unique_acnvs = tool1_acnvs - tool2_acnvs
        tool2_unique_acnvs = tool2_acnvs - tool1_acnvs

        # Save the unique array CNVs
        if len(tool1_unique_acnvs) > 0:
            unique_array_cnvs["tool1"][samplename] = tool1_unique_acnvs
        if len(tool2_unique_acnvs) > 0:
            unique_array_cnvs["tool2"][samplename] = tool2_unique_acnvs
    return unique_array_cnvs


def determine_unique_acnv_types(unique_arraycnvs, arraycnvdata):
    """Determine and return the type of each uniquely found array CNV

    Parameters
    ----------
    unique_arraycnvs : dict
        Uniquely found CNVs for both tools
    arraycnvdata : dict
        Array CNV calls
    """
    unique_acnv_types = {}
    # Calls to `determine_shared_acnv_types` but works for single tools just as well
    tool1_uacnvtypes = determine_shared_acnv_types(unique_arraycnvs["tool1"], arraycnvdata)
    tool2_uacnvtypes = determine_shared_acnv_types(unique_arraycnvs["tool2"], arraycnvdata)
    unique_acnv_types["tool1"] = tool1_uacnvtypes
    unique_acnv_types["tool2"] = tool2_uacnvtypes
    return unique_acnv_types


def compare_fps(tool1_label, tool1_data, tool2_label, tool2_data):
    """Compare the False Positives between classification data of tool1 and tool2.

    Parameters
    ----------
    tool1_label : str
        Name for the first tool
    tool1_data : dict
        Classified CNV calls for the first tool
    tool2_label : str
        Name for the second tool
    tool2_data : dict
        Classified CNV calls for the second tool
    """
    fp_data = {}

    print("[COMPARISON]: Gather False Positive calls for both tools")
    tool1_fps = get_fp_regions(tool1_data)
    tool2_fps = get_fp_regions(tool2_data)

    print("[COMPARISON]: Determine total number of False Positive calls per tool")
    fp_data["Total fps"] = get_total_fps(tool1_fps, tool2_fps)
    print("[COMPARISON]: Determine the number of shared False Positive calls per tool")
    fp_data["Shared fps"] = get_shared_fps(tool1_fps, tool2_fps)
    print("[COMPARISON]: Determine the number of unique False Positive calls per tool")
    fp_data["Unique fps"] = get_unique_fps(tool1_fps, tool2_fps)


def get_fp_regions(tooldata):
    """Get and return False Positive calls for the selected tool.

    Parameters
    ----------
    tooldata : dict
        Tool classification data

    Returns
    -------
    tool_fps : list of str
        False Positives for the tool
    """
    tool_fps = {}
    for samplename in tooldata:
        for cnvcall in tooldata[samplename]:
            if cnvcall.classification== "False Positive":
                if samplename not in tool_fps:
                    tools_fps[samplename] = []
                tool_fps[samplename].append(cnvcall.get_region_str())
    return tool_fps


def get_total_fps(tool1fps, tool2fps):
    fps1 = 0
    for samplename in tool1fps:
        fps1 += len(tool1fps[samplename])

    fps2 = 0
    for samplename in tool2fps:
        fps2 += len(tool2fps[samplename])
    return [fps1, fps2]


def get_shared_fps(tool1data, tool2data):
    shared_fps = {}
    for samplename in set(tool1data.keys()) & set(tool2data.keys()):
        shared_fps[samplename] = set(tool1data[samplename]) & set(tool2data[samplename])
    return shared_fps


def get_unique_fps(tool1data, tool2data):
    unique_fps = {}
    unique_fps["tool1"] = {}
    unique_fps["tool2"] = {}

    for samplename in set(tool1data.keys()) & set(tool2data.keys()):
        unique_fps["tool1"][samplename] = set(tool1data[samplename]) - set(tooldata2[samplename])
        unique_fps["tool2"][samplename] = set(tool2data[samplename]) - set(tooldata1[samplename])
    return unique_fps
