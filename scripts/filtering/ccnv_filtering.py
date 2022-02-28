#!/usr/bin/env python
def filter_classification_results(infileloc, ccnv_data, outfileloc):
    """

    Parameters
    ----------
    infileloc : str
    ccnv_data : dict
    outfileloc : str
        Path to write filtered classifications output file to
    """
    try:
        outfile = open(outfileloc, 'w')
        with open(infileloc, 'r') as infile:
            outfile.write(next(infile))
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                gatkregion = get_gatk_region(filelinedata[1])
                
                if gatkregion[0] in ccnv_data:
                    overlapcount = count_ccnv_overlaps(gatkregion, ccnv_data[gatkregion[0]])
                    if overlapcount < 10:
                        outfile.write(fileline)
        outfile.close()
    except IOError:
        print("Something went wrong :(")


def get_gatk_region(gatkregionstr):
    gatk_chrom = gatkregionstr.split(":")[0]
    gatk_start = int(gatkregionstr.split(":")[1].strip("-")[0])
    gatk_end = int(gatkregionstr.split(":")[1].strip("-")[1])
    return [gatk_chrom, gatk_start, gatk_end]


def count_ccnv_overlaps(gatkregion, ccnvs):
    """Determine the number of overlapping Common CNVs.

    Parameters
    ----------
    gatkregion : str
        GATK4 CNV call region as chr:start-end
    ccnvs : dict
        Common CNVs

    Returns
    -------
    overlap_num : int
        Number of overlapping Common CNVs
    """
    overlap_num = 0
    for ccnv in ccnvs:
        if ccnv_overlaps(gatkregion[1], gatkregion[2], ccnv.cnvstart, ccnv.cnvend):
            overlap_num += 1
    return overlap_num


def ccnv_overlaps(gatkstart, gatkend, ccnvstart, ccnvend):
    """Determine and return whether a specified GATK4 CNV call overlaps with a specified common CNV.

    Parameters
    ----------
    gatkstart : int
        GATK4 CNV call start position
    gatkend : int
        GATK4 CNV call end position
    ccnvstart : int
        Common CNV start position
    ccnvend : int
        Common CNV end position

    Returns
    -------
    bool
        True if GATK4 and common CNV call overlap ; False if not
    """
    return gatkstart <= ccnvend and ccnvstart <= gatkend
