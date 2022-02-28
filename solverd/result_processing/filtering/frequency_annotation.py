#!/usr/bin/env python
import argparse
import statistics
from ccrscall import CcrsCall
from read_combined_ccrs import read_combined_ccrs


def get_params():
    """Define, receive and return set CLI parameter values."""
    freq_annot_args = argparse.ArgumentParser()
    freq_annot_args.add_argument("-i", "--infile", type=str, required=True, dest="infile", help="Path to input file")
    freq_annot_args.add_argument("-n", "--numofsamples", type=int, dest="numofsamples", required=True, help="Total number of processed samples for the batch")
    freq_annot_args.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", help="Path to output directory")
    freq_annot_args.add_argument("-op", "--outprefix", type=str, required=True, dest="outprefix", help="")
    freq_annot_args.add_argument("-p", "--percent-overlap", type=float, dest="percent-overlap", default=80.0, help="")
    freq_annot_args.add_argument("-v", "--verbose", action="store_true", dest="verbose", help="")
    freq_annot_args.add_argument("-c", "--minimum-calls", type=int, dest="minimum-calls", default=3, help="")
    freq_annot_args.add_argument("-s", "--minimum-percentage-samples", type=float, dest="minimum-percentage-samples", default=5.00, help="")
    freq_annot_args.add_argument("-!", "--old-method", action="store_true", dest="old-method", help="Perform the old method of frequency annotation")
    freq_annot_args.add_argument("-f", "--filter-commonxp", action="store_true", dest="filter-commonxp", help="Filter out calls in common xp groups")
    freq_annot_args.add_argument("-m", "--maximum-size", type=int, dest="maximum-size", default=5000000, help="Maximum call size to use to include in groups")
    return vars(freq_annot_args.parse_args())


def determine_call_occurrences(ccrscalls):
    """Determine the occurrence of each specific call.

    Parameters
    ----------
    ccrscalls : dict
        CCRS calls per sample per chromosome.

    Returns
    -------
    call_occurrences : dict
        Occurrences for each call region (as chrom:start-stop)
    """
    call_occurrences = {'+': {}, '-':{}}
    for samplename in ccrscalls:
        for chromname in ccrscalls[samplename]:
            for ccrscall in ccrscalls[samplename][chromname]:
                ccrscalltype = ccrscall.ccrs_call
                if ccrscall.get_region_string() not in call_occurrences[ccrscalltype]:
                    call_occurrences[ccrscalltype][ccrscall.get_region_string()] = 0
                call_occurrences[ccrscalltype][ccrscall.get_region_string()] += 1
    return call_occurrences


def set_call_occurrence_frequency(ccrsdata, calloccurrences, numberofsamples):
    """Set the occurrence and frequency of each call.

    Parameters
    ----------
    ccrsdata : dict
        CCRS calls per sample, per chromosome
    calloccurrences : dict
        Occurrence of each unique call per dup/del per chromosome
    numberofsamples : int
        The number of samples
    """
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                occurr = calloccurrences[ccrscall.ccrs_call][ccrscall.get_region_string()]
                ccrscall.ccrs_occurrence = f"{occurr}/{numberofsamples}"
                ccrscall.ccrs_frequency = round((occurr/numberofsamples)*100, 2)



# =====
# =====METHODS FOR IDENTIFYING CALL GROUPS=====
# =====
def identify_call_groups(ccrscalls, dupdel, maxsize):
    """Identify and return the call groups using the start position.

    Parameters
    ----------
    ccrscalls : dict
        CCRS calls per sample per chromosome
    dupdel : str
        Whether to collect duplications (+) or deletions (-)

    Returns
    -------
    callgroups : dict
        Formed call groups by start position per chromosome
    """
    callgroups = {}
    for samplename in ccrscalls:
        for chromname in ccrscalls[samplename]:
            for ccrscall in ccrscalls[samplename][chromname]:
                if ccrscall.get_call_length() < maxsize:
                    if ccrscall.ccrs_call == dupdel:
                        if chromname not in callgroups:
                            callgroups[chromname] = {}
                        if ccrscall.ccrs_start not in callgroups[chromname]:
                            callgroups[chromname][ccrscall.ccrs_start] = []
                        callgroups[chromname][ccrscall.ccrs_start].append(ccrscall)
    return callgroups


def identify_call_groups_2(ccrscalls, dupdel, maxsize):
    """Identify and return the call groups using the stop position.

    Parameters
    ----------
    ccrscalls : dict
        All CCRS calls per sample per chromosome
    dupdel : str
        Whether to select duplication (+) or deletion (-) calls

    Returns
    -------
    callgroups : dict
        Formed call groups by stop position per chromosome
    """
    callgroups = {}
    for samplename in ccrscalls:
        for chromname in ccrscalls[samplename]:
            for ccrscall in ccrscalls[samplename][chromname]:
                if ccrscall.get_call_length() < maxsize:
                    if ccrscall.ccrs_call == dupdel:
                        if chromname not in callgroups:
                            callgroups[chromname] = {}
                        if ccrscall.ccrs_end not in callgroups[chromname]:
                            callgroups[chromname][ccrscall.ccrs_end] = []
                        callgroups[chromname][ccrscall.ccrs_end].append(ccrscall)
    return callgroups



# =====
# =====METHODS FOR SORTING THE MERGED CALL GROUPS AND CHOOSING THE MEDIAN CALL=====
# =====
def sort_callgroups(callgroups, reversesort):
    """Sort the call groups based on call length."""
    for chromname in callgroups:
        for callstart in callgroups[chromname]:
            callgroups[chromname][callstart].sort(key=lambda x: x.get_call_length(), reverse=reversesort)
    return callgroups


def get_callgroup_representatives(callgroups):
    """Determine and return the call group median."""
    callgroupreps = {}
    for chromname in callgroups:
        if chromname not in callgroupreps:
            callgroupreps[chromname] = {}

        for callstart in callgroups[chromname]:
            if len(callgroups[chromname][callstart]) > 1:
                groupreppos = round(len(callgroups[chromname][callstart])/2)-1
                callgroupreps[chromname][callstart] = callgroups[chromname][callstart][groupreppos]
    return callgroupreps


def form_proper_call_groups(ccrscalls, callgroupreps, minreqoverlap, isstart, isdup):
    """Form the proper call groups, all calls overlapping at least x% with a group representative.

    Parameters
    ----------
    ccrscalls : dict
        All CCRS call per sample, per chromosome
    callgroupreps : dict
        Representatives for callstart or callstop groups
    minreqoverlap : float
        Minimum required overlap to be placed in the group
    """
    grouppref1 = "s" if isstart else "e"
    grouppref2 = "i" if isdup else "d"
    groupnum = 1

    proper_call_groups = {}
    for chromname in callgroupreps:
        for callpos in callgroupreps[chromname]:
            proper_call_groups[f"{grouppref1}{grouppref2}{groupnum}"] = form_proper_single_call_group(ccrscalls, callgroupreps[chromname][callpos], minreqoverlap)
            groupnum +=1
    return proper_call_groups


def form_proper_single_call_group(ccrscalls, callgrouprep, minreqoverlap):
    """Form a proper call group by selecting all calls overlapping with a grouprep."""
    single_group = []
    for samplename in ccrscalls:
        if callgrouprep.ccrs_chrom in ccrscalls[samplename]:
            for ccrscall in ccrscalls[samplename][callgrouprep.ccrs_chrom]:
                if ccrscall.ccrs_call == callgrouprep.ccrs_call:
                    if determine_call_rep_overlap(ccrscall.ccrs_start, ccrscall.ccrs_end, callgrouprep.ccrs_start, callgrouprep.ccrs_end) >= minreqoverlap:
                        single_group.append(ccrscall)
    return single_group


def determine_call_rep_overlap(callstart, callend, repstart, repend):
    """Determine the percentage overlap between a group call and the group representative.

    Parameters
    ----------
    callstart : int
        Starting position of the call
    callend : int
        Ending position of the call
    repstart : int
        Starting position of the representative call
    repend : int
        Ending position of the representative call
    """
    call_length = callend - callstart
    overlapnum = max(0, min(repend, callend) - max(repstart, callstart))
    return round((overlapnum/call_length)*100, 2)


def calculate_xpgroup_occurrences(xpgroups, calloccurrences):
    """Calculate and return the group occurrences.

    Parameters
    ----------
    xpgroups : dict
        CCRS calls in groups with x% overlap
    calloccurrences : dict
        Occurrence for each unqiue call

    Returns
    -------
    group_occurrences : dict
        Total occurrence for each xp group
    """
    group_occurrences = {}
    for groupname in xpgroups:
        group_occurrence = 0
        calls_processed = []
        for ccrscall in xpgroups[groupname]:
            if ccrscall.get_region_string() not in calls_processed:
                group_occurrence += calloccurrences[ccrscall.get_region_string()]
                calls_processed.append(ccrscall.get_region_string())
        group_occurrences[groupname] = group_occurrence
    return group_occurrences


def determine_common_xpgroups(xpgroupoccurrences, numofsamples, minoccurrence, minfrequency):
    """Determine whether the xp groups are common or not based on the occurrence and frequency.

    Parameters
    ----------
    xpgroupoccurrences : dict
        Occurrence for each xp group
    numofsamples : int
        Number of samples for the BED (batch)
    minoccurrence : int
        Minimum occurrence required (default: 3)
    minfrequency : float
        Minimum frequency required (default: 5.0)
    """
    xpcommons = {}
    for xpname in xpgroupoccurrences:
        xpcommons[xpname] = False
        xpfrequency = round((xpgroupoccurrences[xpname]/numofsamples)*100, 2)
        if xpgroupoccurrences[xpname] >= minoccurrence and xpfrequency >= minfrequency:
            xpcommons[xpname] = True
    return xpcommons


def callgroups_is_common(callgroup, totalsamples, minnumofcalls, minpercofsamples):
    """Determine whether the calls in the call.

    Parameters
    ----------
    callgroup : dict
    totalsamples : int
    minnumofcalls : int
        
    minpercofsamples : float
        Minimum percentage of samples that a call group needs
    """
    if len(callgroup) >= minnumofcalls:
        numofgroupsamples = len(set([x.ccrs_sample for x in callgroup]))
        groupsampleperc = round( (numofgroupsamples/totalsamples)*100, 2)
        if groupsampleperc >= minpercofsamples:
            # Set the ccrscall as common
            for ccrscall in callgroup:
                ccrscall.callgroup_common = True



# =====
# =====METHODS FOR SETTING THE GROUP OCCURRENCES AND FREQUENCIES=====
# =====

def set_callgroups_occurrences_frequencies(callgroups, calloccurrences, numberofsamples):
    """Determine the group occurrences and frequencies for all callstart/callstop groups.

    Parameters
    ----------
    callgroups : dict
        Call start/stop groups of dup/del calls to determine occurrences and frequencies for
    calloccurrences : dict
        Call occurrences for dup/del
    numberofsamples : int
        The processed number of samples
    """
    for chromname in callgroups:
        for callposition in callgroups[chromname]:
            group_occur_freq = get_callgroup_occurrence_frequency(callgroups[chromname][callposition], calloccurrences, numberofsamples)
            for ccrscall in callgroups[chromname][callposition]:
                ccrscall.ccrs_callgroup_occurrence = group_occur_freq[0]
                ccrscall.ccrs_callgroup_frequency = group_occur_freq[1]


def get_callgroup_occurrence_frequency(callgroup, calloccurrences, numberofsamples):
    """Determine the occurrence and frequency of a single callstart/callstop group.

    Parameters
    ----------
    callgroup : list of CcrsCall
        CCRS calls of a single group
    calloccurrences : dict
        Occurrences of calls
    numberofsamples : int
        Number of processed samples
    """
    occurr = 0
    for ccrscall in callgroup:
        if ccrscall.get_region_string() in calloccurrences:
            occurr += calloccurrences[ccrscall.get_region_string()]
    return [f"{occurr}/{numberofsamples}", round((occurr/numberofsamples)*100, 2)]


def set_xpgroup_occurrences_frequencies(xpgroups, xpoccurrences, numofsamples):
    """Set the highest group occurrence and frequency for each CCRS call.

    Parameters
    ----------
    xpgroups : dict
        The xp groups
    xpoccurrences : dict
        Occurrences for each xp group
    numofsamples : int
        Total number of samples for the batch
    """
    for groupname in xpgroups:
        group_freq = round((xpoccurrences[groupname]/numofsamples)*100, 2)
        group_occurr = xpoccurrences[groupname]
        for ccrscall in xpgroups[groupname]:
            if group_freq > ccrscall.ccrs_callgroup_frequency:
                ccrscall.ccrs_callgroup_occurrence = f"{group_occurr}/{numofsamples}"
                ccrscall.ccrs_callgroup_frequency = group_freq
                ccrscall.ccrs_callgroup_name = groupname


def set_non_xpgroup_occurrences_frequencies(ccrscalls):
    """Set the group occurrences and frequencies for calls that were not part of an xp group."""
    groupnum = 1
    for samplename in ccrscalls:
        for chromname in ccrscalls[samplename]:
            for ccrscall in ccrscalls[samplename][chromname]:
                if ccrscall.ccrs_callgroup_occurrence == "":
                    ccrscall.ccrs_callgroup_name = f"xi{groupnum}" if ccrscall.ccrs_call == '+' else f"xd{groupnum}"
                    ccrscall.ccrs_callgroup_occurrence = ccrscall.ccrs_occurrence
                    ccrscall.ccrs_callgroup_frequency = ccrscall.ccrs_frequency
                    groupnum += 1


def determine_common_xpgroup_calls(callgroups, minnumofcalls, minpercofsamples):
    """Determine which calls can be considered common by ."""



# =====
# =====DO THE MAIN STUFF=====
# =====
def main():
    freq_annot_params = get_params()
    freq_annot_params["outdir"] = freq_annot_params["outdir"]+"/" if not freq_annot_params["outdir"].endswith("/") else freq_annot_params["outdir"]

    # Obtain the CCRS data.
    print("[-READING THE CCRS DATA-]")
    ccrs_data = read_combined_ccrs(freq_annot_params["infile"])
    ccrs_header = ccrs_data[0]
    ccrs_calls = ccrs_data[1]

    # Determine the call occurrences
    print("[-DETERMINING CALL OCCURRENCES-]")
    call_occurrences = determine_call_occurrences(ccrs_calls)

    if freq_annot_params["verbose"]:
        print("...Call occurrences...")
        tmp_show_occurrences(call_occurrences)
        print("")


    # Set the call occurrence and frequency
    print("[-SETTING THE CALL OCCURRENCES-]")
    set_call_occurrence_frequency(ccrs_calls, call_occurrences, freq_annot_params["numofsamples"])


    # Determine the call groups based on the same start position
    print("[-DETERMINING CALL GROUPS BY START AND END POSITION-]")
    dup_start_call_groups = identify_call_groups(ccrs_calls, '+', freq_annot_params["maximum-size"])
    del_start_call_groups = identify_call_groups(ccrs_calls, '-', freq_annot_params["maximum-size"])
    dup_stop_call_groups = identify_call_groups_2(ccrs_calls, '+', freq_annot_params["maximum-size"])
    del_stop_call_groups = identify_call_groups_2(ccrs_calls, '-', freq_annot_params["maximum-size"])

    if freq_annot_params["verbose"]:
        print("...Duplication callstart groups...")
        tmp_show_callgroups(dup_start_call_groups)
        print("")
        print("...Deletion callstart groups...")
        tmp_show_callgroups(del_start_call_groups)
        print("")
        print("...Duplication callstop groups...")
        tmp_show_callgroups(dup_stop_call_groups)
        print("")
        print("...Deletion callstop groups...")
        tmp_show_callgroups(del_stop_call_groups)
        print("")


    # Sort the calls based on length (smallest -> longest)
    print("[-SORTING DUPLICATION AND DELETION CALL GROUPS-]")
    dup_start_call_groups = sort_callgroups(dup_start_call_groups, False)
    del_start_call_groups = sort_callgroups(del_start_call_groups, False)
    dup_stop_call_groups = sort_callgroups(dup_stop_call_groups, False)
    del_stop_call_groups = sort_callgroups(del_stop_call_groups, False)

    if freq_annot_params["verbose"]:
        print("...Sorted duplication callstart groups...")
        tmp_show_sorted_callgroups(dup_start_call_groups)
        print("")
        print("...Sorted deletion callstart groups")
        tmp_show_sorted_callgroups(del_start_call_groups)
        print("")
        print("...Sorted duplication callstop groups...")
        tmp_show_sorted_callgroups(dup_stop_call_groups)
        print("")
        print("...Sorted deletion callstop groups...")
        tmp_show_sorted_callgroups(del_stop_call_groups)
        print("")


    # Determine the callstart and callstop group representatives (by selecting the median call based on size)
    print("[-DETERMINING CALL GROUP REPRESENTATIVES-]")
    dup_callstart_representatives = get_callgroup_representatives(dup_start_call_groups)
    dup_callstop_representatives = get_callgroup_representatives(dup_stop_call_groups)
    del_callstart_representatives = get_callgroup_representatives(del_start_call_groups)
    del_callstop_representatives = get_callgroup_representatives(del_stop_call_groups)

    if freq_annot_params["verbose"]:
        print("...Duplication callstart group representatives...")
        tmp_show_group_representatives(dup_callstart_representatives)
        print("")
        print("...Deletion callstart group representatives...")
        tmp_show_group_representatives(del_callstart_representatives)
        print("")
        print("...Duplication callstop group representatives...")
        tmp_show_group_representatives(dup_callstop_representatives)
        print("")
        print("...Deletion callstop group representatives...")
        tmp_show_group_representatives(del_callstop_representatives)
        print("")


    # Determine the
    print("[-DETERMINING THE PROPER DUPLICATION AND DELETION CALL GROUPS BASED ON OVERLAP WITH THE GROUP REPRESENTATIVES-]")
    xp_dup_callstart = form_proper_call_groups(ccrs_calls, dup_callstart_representatives, freq_annot_params["percent-overlap"], True, True)
    xp_del_callstart = form_proper_call_groups(ccrs_calls, del_callstart_representatives, freq_annot_params["percent-overlap"], True, False)
    xp_dup_callstop = form_proper_call_groups(ccrs_calls, dup_callstop_representatives, freq_annot_params["percent-overlap"], False, True)
    xp_del_callstop = form_proper_call_groups(ccrs_calls, del_callstop_representatives, freq_annot_params["percent-overlap"], False, False)


    # Calculate the xp group occurrences
    print("[-CALCULATE THE XP GROUP OCCURRENCES-]")
    xp_dupstart_occurrences = calculate_xpgroup_occurrences(xp_dup_callstart, call_occurrences['+'])
    xp_delstart_occurrences = calculate_xpgroup_occurrences(xp_del_callstart, call_occurrences['-'])
    xp_dupstop_occurrences = calculate_xpgroup_occurrences(xp_dup_callstop, call_occurrences['+'])
    xp_delstop_occurrences = calculate_xpgroup_occurrences(xp_del_callstop, call_occurrences['-'])

    if freq_annot_params["verbose"]:
        print("...Duplication...")
        tmp_show_xpgroup_occurrences(xp_dupstart_occurrences)
        tmp_show_xpgroup_occurrences(xp_delstart_occurrences)
        tmp_show_xpgroup_occurrences(xp_dupstop_occurrences)
        tmp_show_xpgroup_occurrences(xp_delstop_occurrences)


    # Set the call group occurrences, frequencies and names
    print("[-SET THE GROUP OCCURRENCE AND FREQUENCY FOR CALLS IN XP GROUPS-]")
    set_xpgroup_occurrences_frequencies(xp_dup_callstart, xp_dupstart_occurrences, freq_annot_params["numofsamples"])
    set_xpgroup_occurrences_frequencies(xp_del_callstart, xp_delstart_occurrences, freq_annot_params["numofsamples"])
    set_xpgroup_occurrences_frequencies(xp_dup_callstop, xp_dupstop_occurrences, freq_annot_params["numofsamples"])
    set_xpgroup_occurrences_frequencies(xp_del_callstop, xp_delstop_occurrences, freq_annot_params["numofsamples"])
    # tmp_show_ccrs_group_occurr(ccrs_calls)


    # Set the 
    print("[-SETTING THE GROUP OCCURRENCE AND FREQUENCY FOR CALLS NOT IN ANY XP GROUP-]")
    set_non_xpgroup_occurrences_frequencies(ccrs_calls)
    tmp_show_ccrs_group_occurr(ccrs_calls)


    # Determine which XP groups are common according to the set minimum-calls and minimum-percentage-samples
    print("[-DETERMINING WHICH XP GROUPS ARE COMMON-]")
    common_xp_groups = {}
    common_xp_groups.update(determine_common_xpgroups(xp_dupstart_occurrences, freq_annot_params["numofsamples"], freq_annot_params["minimum-calls"], freq_annot_params["minimum-percentage-samples"]))
    common_xp_groups.update(determine_common_xpgroups(xp_delstart_occurrences, freq_annot_params["numofsamples"], freq_annot_params["minimum-calls"], freq_annot_params["minimum-percentage-samples"]))
    common_xp_groups.update(determine_common_xpgroups(xp_dupstop_occurrences, freq_annot_params["numofsamples"], freq_annot_params["minimum-calls"], freq_annot_params["minimum-percentage-samples"]))
    common_xp_groups.update(determine_common_xpgroups(xp_delstop_occurrences, freq_annot_params["numofsamples"], freq_annot_params["minimum-calls"], freq_annot_params["minimum-percentage-samples"]))


    # Write the new combined CCRS file
    print("[-WRITING CCRS CALLS TO FILE-]")
    outfilepath = freq_annot_params["outdir"] + freq_annot_params["outprefix"] + ".called.seg"
    wrote_file = write_ccrs_calls(outfilepath, ccrs_calls, ccrs_header, common_xp_groups, freq_annot_params["filter-commonxp"])
    print(f"...Wrote new combined CCRS file?: {wrote_file}...")


    # Write the group data to file(s)
    print("[-WRITING THE XP GROUP CONTENTS TO FILE-]")
    dup_outpref = freq_annot_params["outdir"] + freq_annot_params["outprefix"] + "_dup_"
    del_outpref = freq_annot_params["outdir"] + freq_annot_params["outprefix"] + "_del_"
    write_xp_groups(f"{dup_outpref}callstart.txt", xp_dup_callstart)
    write_xp_groups(f"{del_outpref}callstart.txt", xp_del_callstart)
    write_xp_groups(f"{dup_outpref}callstop.txt", xp_dup_callstop)
    write_xp_groups(f"{del_outpref}callstop.txt", xp_del_callstop)



# =====TEMPORARY METHODS TO DISPLAY DATA=====
def tmp_show_occurrences(occurrencedata):
    """Display the call occurrences."""
    for dupdel in occurrencedata:
        print(f"[{dupdel}]")
        for callregion in occurrencedata[dupdel]:
            print(f"{callregion}\t{occurrencedata[dupdel][callregion]}")


def tmp_show_callgroups(callgroups):
    """Display the formed call groups."""
    for chromname in callgroups:
        for callstart in callgroups[chromname]:
            aap = f"{chromname}:{callstart}:\t"
            for ccrscall in callgroups[chromname][callstart]:
                aap += ccrscall.get_region_string()
            print(aap)


def tmp_show_group_representatives(callgroupreps):
    """Display the group representatives."""
    for chromname in callgroupreps:
        for callstart in callgroupreps[chromname]:
            print(callgroupreps[chromname][callstart].get_region_string())


def tmp_show_sorted_callgroups(callgroups):
    for chromname in callgroups:
        for callstart in callgroups[chromname]:
            for ccrscall in callgroups[chromname][callstart]:
                print(f"{ccrscall.ccrs_chrom}:{ccrscall.ccrs_start}-{ccrscall.ccrs_end}")


def tmp_show_xp_groups(xpgroups):
    for groupname in xpgroups:
        for ccrscall in xpgroups[groupname]:
            print(f"{groupname}\t({ccrscall.ccrs_sample}){ccrscall.get_region_string()}\n")


def tmp_show_xpgroup_occurrences(xpgroupsoccurrs):
    for groupname in xpgroupsoccurrs:
        print(f"{groupname}: {xpgroupsoccurrs[groupname]}")


def tmp_show_ccrs_group_occurr(ccrscalls):
    for samplename in ccrscalls:
        for chromname in ccrscalls[samplename]:
            for ccrscall in ccrscalls[samplename][chromname]:
                print(f"{ccrscall.ccrs_callgroup_name}: {ccrscall.ccrs_callgroup_occurrence}\t{ccrscall.ccrs_call}\t{ccrscall.ccrs_occurrence}\t{ccrscall.get_region_string()}\t{ccrscall.get_call_length()}")


def tmp_show_occurrs(ccrscalls):
    for samplename in ccrscalls:
        for chromname in ccrscalls[samplename]:
            for ccrscall in ccrscalls[samplename][chromname]:
                print(ccrscall.ccrs_occurrence)


def write_ccrs_calls(outfileloc, ccrscalls, ccrsheader, xpcommons, filtercommons):
    """Write the CCRS calls with their individual and xp group frequencies.

    Parameters
    ----------
    outfileloc : str
        Path to write output file to
    ccrscalls : dict
        CCRS calls per sample per chromosome
    ccrsheader : str
        Original combined CCRS header to expand for the new file
    xpcommons : dict
        Which xp groups are and are not common
    filtercommons : bool
        Whether to filter (not write) calls in common xp groups
    """
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(f"{ccrsheader.strip()}\tCall_Occurrence\tCall_Frequency\tCall_Group\tGroup_Occurrence\tGroup_Frequency\n")
            for samplename in ccrscalls:
                for chromname in ccrscalls[samplename]:
                    for ccrscall in ccrscalls[samplename][chromname]:
                        if filtercommons:
                            if ccrscall.ccrs_callgroup_name in xpcommons:
                                if not xpcommons[ccrscall.ccrs_callgroup_name]:
                                    outfile.write(ccrscall.to_ccrs_file_line_2())
                            else:
                                outfile.write(ccrscall.to_ccrs_file_line_2())
                        else:
                            outfile.write(ccrscall.to_ccrs_file_line_2())
        file_written = True
    except IOError:
        print("Could not write new combined CCRS file")
    finally:
        return file_written


def write_xp_groups(outfileloc, xpgroups):
    """Write proper percentage overlap groups to file.

    Parameters
    ----------
    outfileloc : str
        Path to write output file to
    xpgroups : dict
        Duplication/Deletion groupcalls to write to file

    Return
    ------
    file_written : bool
        True if file has been successfully written, False if not
    """
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Groupname\tSample_Call\n")
            for groupname in xpgroups:
                for ccrscall in xpgroups[groupname]:
                    outfile.write(f"{groupname}\t({ccrscall.ccrs_sample}){ccrscall.get_region_string()}\n")
            file_written = True
    except IOError:
        print("Could not write output file")
    finally:
        return file_written


if __name__ == "__main__":
    main()
