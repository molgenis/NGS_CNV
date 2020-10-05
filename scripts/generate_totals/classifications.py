def generate_classification_totals(classificationdata):
    classificationtotals = {}
    for samplename in classificationdata:
        for gatkcall in classificationdata[samplename]:
            label = determine_clasification_totals_label(gatkcall.classification)
            if label not in classificationtotals:
                classificationtotals[label] = 0
            classificationtotals[label] = classificationtotals[label] + 1
    return classificationtotals


def determine_clasification_totals_label(classification):
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
    return classlabel
