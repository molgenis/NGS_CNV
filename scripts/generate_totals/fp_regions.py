from classes.fpoverlap import FpOverlap
import shared_methods.shared_methods as smsm

def get_duplicated_regions(fpregioncounts):
    """Fetch and return false positive regions

    Parameters
    ----------
    fpregioncounts : dict
        Occurence counts per False Positive region

    Returns
    -------
    dup_fpregions : list of str
        List of false postive regions with count > 1
    """
    dup_fpregions = {}
    for fpregion in fpregioncounts:
        if fpregioncounts[fpregion] > 1:
            dup_fpregions[fpregion] = fpregioncounts[fpregion]
    return dup_fpregions


def get_unique_regions(fpregioncounts):
    """Fetch and return unique regions based on the counts.

    Parameters
    ----------
    fpregioncounts : dict
        Number of occurences for each region

    Return
    ------
    uni_fpregions : list of str
        List of unique false positive regions with count == 1
    """
    uni_fpregions = {}
    for fpregion in fpregioncounts:
        if fpregioncounts[fpregion] == 1:
            uni_fpregions[fpregion] = fpregioncounts[fpregion]
    return uni_fpregions


def determine_similar_regions(dupregions, uniregions, min_req_overlap):
    """Determine which False Positive regions have at least a minimum required overlap.

    Parameters
    ----------
    dupregions : dict
        
    uniregions : dict
    min_req_overlap : int
        Minimum percentage overlap

    Returns
    -------
    simfpregions : dict
        False Positive regions with 
    """
    simfpregions = {}
    simfpregions = __similar_with_duplicates__(dupregions, uniregions, simfpregions, min_req_overlap)
    simfpregions = __similar_with_uniques__(uniregions, simfpregions, min_req_overlap)
    return simfpregions


def __similar_with_duplicates__(dupregions, uniregions, simregions, minpercoverlap):
    """Determine similar 

    Parameters
    ----------
    dupregions : dict
    uniregions : dict
    simregions : dict
    minpercoverlap : int
    """
    for uregion in uniregions:
        for dregion in dupregions:
            fpoverlap = regions_overlap(uregion, dregion, simregions, minpercoverlap)

            if fpoverlap is not None:
                if uregion not in simregions:
                    simregions[uregion] = FpOverlap(uregion)
                simregions[uregion].add_overlap(fpoverlap[0], fpoverlap[1])
    return simregions


def __similar_with_uniques__(uniregions, simregions, minpercoverlap):
    for uregion in uniregions:
        for ouregion in uniregions:
            fpoverlap = regions_overlap(uregion, ouregion, simregions, minpercoverlap)

            if fpoverlap is not None:
                if uregion not in simregions:
                    simregions[uregion] = FpOverlap(uregion)
                simregions[uregion].add_overlap(fpoverlap[0], fpoverlap[1])
    return simregions


def regions_overlap(selected_region, other_region, simregions, min_overlap_percentage):
    """Check and return whether a selected_region overlaps with the other_region for at least the set overlap_percentage.

    Parameters
    ----------
    selected_region : str
    other_region : str
    min_overlap_percentage : int
        Minimum required percentage overlap

    Returns
    -------
    list of str and int or None
        Overlapping region and percentage overlap ; None if there is no overlap
    """
    selected_region_data = smsm.get_gatk_region(selected_region)
    other_region_data = smsm.get_gatk_region(other_region)

    if selected_region != other_region:
        if selected_region_data[0] == other_region_data[0]:
            if selected_region_data[1] <= other_region_data[2] and other_region_data[1] <= selected_region_data[2]:
                print(f"{selected_region} :: {other_region}")
                print(f"{selected_region_data} :: {other_region_data}")

                if not similarfp_already_found(selected_region, other_region, simregions):
                    selected_region_size = selected_region_data[2] - selected_region_data[1]

                    # Determine the size of the overlap
                    non_overlap_size = 0
                    if selected_region_data[1] < other_region_data[1]:
                        non_overlap_size += other_region_data[1] - selected_region_data[1]
                    if selected_region_data[1] > other_region_data[2]:
                        non_overlap_size += selected_region_data[2] - other_region_data[2]

                    # Determine the overlap percentage
                    overlap_size = selected_region_size - non_overlap_size
                    if overlap_size > 0:
                        perc_overlap = (overlap_size/selected_region_size)*100
                        if perc_overlap >= min_overlap_percentage:
                            return [other_region, perc_overlap]
    return None


def similarfp_already_found(selectedregion, otherregion, simfpregions):
    """Check and return whether a similar fpregion has already been found.

    This is to avoid an A<->B entry to be recorded as B<->A as well.

    Parameters
    ----------
    selectedregion : str
    otherregion : str
    """
    if otherregion in simfpregions:
        return selectedregion in simfpregions[otherregion].overlaps
    return False


def filter_uniquefps_with_similars(simfp_filterlist, unifpregions):
    """Filter similar False Positive regions from the set of unique False Positive regions.

    Parameters
    ----------
    simfpregions : list
        Similar False Positive regions
    unifpregions : dict
        Unique False Positive regions
    """
    for simregion in simfp_filterlist:
        if simregion in unifpregions:
            del unifpregions[simregion]
    return unifpregions


def get_similar_filterlist(simfpregions):
    """Construct and return the list of all similar False Positive regions.

    Parameters
    ----------
    simfpregions : dict
        Similar False Positive regions
    """
    simfilterlist = []
    for simregion in simfpregions:
        simfilterlist.append(simregion)
        simfilterlist.extend(list(simfpregions[simregion].overlaps.keys()))
    return list(set(simfilterlist))
