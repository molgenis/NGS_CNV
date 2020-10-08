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
    """Fetch and return 

    Parameters
    ----------
    fpregioncounts : dict
        

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
    simfpregions = {}
    simfpregions = self.__similar_with_duplicates__(dupregions, uniregions, simfpregions, min_req_overlap)
    simfpregions = self.__similar_with_uniques__(uniregions, simfpregions, min_req_overlap)
    return simfpregions


def __similar_with_duplicates__(dupregions, uniregions, simregions, minpercoverlap):
    for uregion in uniregions:
        for dregion in dupregions:
            fpoverlap = regions_overlap(uregion, dregion, minpercoverlap)

            if fpoverlap is not None:
                if uregion not in simregions:
                    simregions[uregion] = FpOverlap(uregion)
                simregions[uregion].add_overlap(fpoverlap[0], fpoverlap[1])
    return simregions


def __similar_with_uniques__(unireqions, simregions, minpercoverlap):
    for uregion in uniregions:
        for ouregion in uniregions:
            fpoverlap = regions_overlap(uregion, ouregion, minpercoverlap)

            if fpoverlap is not None:
                if uregion not in simregions:
                    simregions[uregion] = FpOverlap(uregion)
                simregions[uregions].add_overlap(fpoverlap[0], fpoverlap[1])
    return simregions


def regions_overlap(selected_region, other_region, min_overlap_percentage):
    """Check and return whether a selected_region overlaps with the other_region for at least the set overlap_percentage."""
    if selected_region[1] <= other_region[2] and other_region[1] <= selected_region[2]:
        selected_region_data = smsm.get_gatk_region(selected_region)
        other_region_data = smsm.get_gatk_region(other_region)
        selected_region_size = selected_region_data[2] - selected_region_data[1]

        # Determine the size of the overlap
        non_overlap_size = 0
        if selected_region_data[1] < other_region_data[1]:
            non_overlap_size += other_region_data[1] - selected_region_data[1]
        if selected_region_data[1] > other_region_data[2]:
            non_overlap_size += selected_region_data[2] - other_region_data[2]

        # Determine the overlap percentage
        overlap_size = selected_region_size - non_overlap_size
        perc_overlap = (overlap_size/selected_region_size)*100

        if perc_overlap >= min_overlap_percentage:
            return [other_region, perc_overlap]
    return None
