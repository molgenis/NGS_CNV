#!/usr/bin/env python
def read_conifer_data(coniferfileloc, sampledata, exondata, probedata):
    conifer_data = {}
    try:
        with open(coniferfileloc, 'r') as coniferfile:
            next(coniferfile)
            for fileline in coniferfile:
                filelinedata = fileline.strip().split("\t")
                cfcall = ConiferCall(filellinedata[0], filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), filelinedata[4])
                cfcall.exons = None
                cfcall.probes = None
                conifer_data
    except IOError:
        print(f"")
    finally:
        return conifer_data


def conifer_evaluate(sample_data, array_cnvs, conifer_cnv_calls):
    """Evaluate the Conifer calls against the array CNV calls.

    Parameters
    ----------
    sample_data : dict
        Sample translation table
    array_cnvs : dict
        Array called CNVs
    conifer_cnv_calls : dict
        Conifer CNV calls per sample
    """
    for samplename in conifer_cnv_calls:
        print(f"Evaluating sample {samplename}")
        for conifercnv in conifer_cnv_calls[samplename]:
            if conifercnv.cnv_sample_pseudo in array_cnvs:
                arraysamplecnvs = arraycnv[conifercnv.cnv_sample_pseudo]
                conifer_classify_cnv(conifercnv, arraysamplecnvs, num_of_exons, num_of_probes)


def conifer_classify_cnv(conifercnv, array_cnvs, num_of_exons, num_of_probes, percent_overlap):
    """Classify a Conifer CNV using the array CNVs.

    Parameters
    ----------
    conifercnv : ConiferCall
        Single Conifer CNV call
    array_cnvs : dict
        Array CNV calls per sample
    num_of_exons : int
        Minimum required number of exons for a CNV call to be considered WES Informative
    num_of_probes: int
        Minimum required number of probes for a CNV call to be considered Array Informative
    percent_overlap : int
        Minimum required percent overlap for a CNV to be considered a True Positive
    """
    match_with_array = False
    classification = ""
    callresult = ""
    for arraycnv in array_cnvs:
        

