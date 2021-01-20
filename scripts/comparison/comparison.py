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
    tool1_classification_numbers = gtc.generate_classification_totals(tool1_data)
    tool2_classification_numbers = gtc.generate_classification_totals(tool2_data)
    classification_numbers = {"True Positive": [tool1_classification_numbers["True Positive"], tool2_classification_numbers["True Positive"]],
                              "False Positive": [tool1_classification_numbers["False Positive"], tool2_classification_numbers["False Positive"]],
                              "Array Non-Informative": [tool1_classification_numbers["Array Non-Informative"], tool2_classification_numbers["Array Non-Informative"]],
                              "WES Non-Informative": [tool1_classification_numbers["WES Non-Informative"], tool2_classification_numbers["WES Non-Informative"]],
                              "Array & WES Non-Informative": [tool1_classification_numbers["Array & WES Non-Informative"], tool2_classification_numbers["Array & WES Non-Informative"]]}
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
    for samplename in shared_arraycnv:
        for shared_acnv in shared_arraycnvs[samplename]:
            if shared_acnv in arraydata[samplename]:
                array_type = arraydata[samplename][shared_acnv].cnv_class
                if array_type not in shared_acnv_types:
                    shared_acnv_types[array_type] = 0
                shared_acnv_types[array_type] += 1
    return shared_acnv_types


def determine_unique_array_cnvs(found_array_cnv_data):
    """Determine which array cNVs are unique for each tool.

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
