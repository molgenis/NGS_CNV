def add_probes(cnv, probes):
    """Identify probes overlapping with the the CNV call.

    Parameters
    ----------
    cnv : Cnv or ConiferCall or ExomeDepthCall
        GATK4 CNV
    probes : dict
        Probe data
    """
    probes_to_add = []
    if cnv.cnv_chrom in probes:
        chromprobes = probes[cnv.cnv_chrom]
        for chromprobe in chromprobes:
            if cnv.segment_overlap(chromprobe.probe_start, chromprobe.probe_end):
                probes_to_add.append(chromprobe)
    return probes_to_add


def add_exons(cnv, exondata):
    """Identify exons overlapping with the CNV call.

    Parameters
    ----------
    cnv : Cnv or ConiferCall or ExomeDepthCall
        GATK4, Conifer or ExomeDepth CNV to 
    exondata : dict
        Exon data containing the exons to add to the CNV call
    """
    exons_to_add = []
    chromexons = exondata[cnv.cnv_chrom]
    for chromexon in chromexons:
        if cnv.segment_overlap(chromexon.exon_start, chromexon.exon_end):
            exons_to_add.append(chromexon)
    return exons_to_add


def evaluate_cnv_calls(cnv_calls, array_cnvs, exon_data, probe_data):
    """Evaluate GATK4, Conifer or ExomeDeph CNV calls with array CNV calls.

    Parameters
    ----------
    cnv_calls
        GATK4, Conifer or ExomeDepth CNV calls to be evaluated
    array_cnvs : dict
        Array CNV calls saved per sample
    exon_data : dict
        Exon data saved per chromosome
    probe_data : dict
        Probe data saved per chromosome

    Returns
    -------
    cnv_calls : dict
        Classified GATK4, Conifer or ExomeDepth CNV calls
    """
    for samplepseudo in cnv_calls:
        for cnvcall in cnv_calls[samplepseudo]:
            if cnvcall.cnv_sample_pseudo in array_cnvs:
                arraysamplecnvs = array_cnvs[cnvcall.cnv_sample_pseudo]
                classify_cnv(cnvcall, arraysamplecnvs)
    return cnv_calls
