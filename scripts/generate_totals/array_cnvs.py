def determine_arraycnvs_found(infileloc, arraycnvs):
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


def determine_arraycnvs_missed(foundcnvs, arraycnvs):
    """Determine and return which array CNVs were missed.

    Parameters
    ----------
    foundcnvs : dict
        Array CNVs that have been found
    arraycnvs : dict
        
    """
    missedcnvs = {}
    for samplename in arraycnvs:
        if samplename not in foundcnvs:
            missedcnvs[samplename] = []
            for acnvregion in arraycnvs[samplename]:
                missedcnvs[samplename].append(arraycnvs[samplename][acnvregion])
        else:
            for acnvregion in arraycnvs[samplename]:
                if not is_arraycnv_in_foundcnvs(acnvregion, foundcnvs[samplename]):
                    if samplename not in missedcnvs:
                        missedcnvs[samplename] = []
                    missedcnvs[samplename].append(arraycnvs[samplename][acnvregion])
    return missedcnvs


def is_arraycnv_in_foundcnvs(arraycnvregion, foundsampledata):
    """Determine and return whether an array CNV region is in the list of already found array CNVs.

    Parameters
    ----------
    arraycnvregion : str
        Array CNV region
    foundsampledata : dict

    Returns
    -------
    acnv_present : bool
        Whether a specific array CNV is present in the set of found array CNV
    """
    acnv_present = False
    for foundcnv in foundsampledata:
        if arraycnvregion == foundcnv.get_region():
            acnv_present = True
    return acnv_present


def summarize_arraycnv_types(arraycnvsfound):
    """Summarize the number of found array CNV types.

    Parameters
    ----------
    arraycnvsfound : dict
        Array CNVs that were found
    """
    typecounts = {}
    for samplename in arraycnvsfound:
        for facnv in arraycnvsfound[samplename]:
            if facnv.cnv_class not in typecounts:
                typecounts[facnv.cnv_class] = 0
            typecounts[facnv.cnv_class] = typecounts[facnv.cnv_class] + 1
    return typecounts
