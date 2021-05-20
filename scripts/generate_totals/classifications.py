def generate_classification_totals(classificationdata, tpperacnv):
    """Determine and return the totals for TP, FP, ANI, WNI, AWNI.

    Parameters
    ----------
    classificationdata : dict
        Classified CNV calls
    tpperacnv : bool
        Count TPs per array CNV or not

    Returns
    -------
    classificationtotals : dict
        Totals for each classification label
    """
    classificationtotals = {}
    for samplename in classificationdata:
        for gatkcall in classificationdata[samplename]:
            label = determine_clasification_totals_label(gatkcall.classification)
            if label not in classificationtotals:
                classificationtotals[label] = 0
            classificationtotals[label] = classificationtotals[label] + 1
    if tpperacnv:
        numoftps = determine_tps_as_arraycnvs(classificationdata)
        classificationtotals["True Positive"] = numoftps
    return classificationtotals


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


def determine_tps_as_arraycnvs(cnvcalls):
    """Only count True Positives as array CNVs.
    Two True Positives overlapping with the same array CNV will be counted as 1 True Positive.

    Parameters
    ----------
    cnvcalls : dict
        Classified CNV calls

    Returns
    -------
    num_of_tps : int
        Number of True Positives
    """
    arraycnvs_found = {}
    num_of_tps = 0
    for samplename in cnvcalls:
        for cnvcall in cnvcalls[samplename]:
            if cnvcall.classification == "POSITIVE" or cnvcall.classification == "TRUE POSITIVE":
                if samplename not in arraycnvs_found:
                    arraycnvs_found[samplename] = []
                if cnvcall.arraycnv not in arraycnvs_found[samplename]:
                    num_of_tps += 1
                    arraycnvs_found[samplename].append(cnvcall.arraycnv)
    return num_of_tps
