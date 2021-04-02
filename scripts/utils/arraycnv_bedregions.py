def get_arraycnv_bedregions(bedfiledata, arraycnvregion):
    """Select and return BED file regions overlapping with an array CNV region.

    Parameters
    ----------
    bedfiledata : dict
        Read BED file regions
    arraycnvregion : str
        Array CNV as a region (chr:start-end)
    """
    acnv_data = arraycnvregion.split(":")
    acnv_chrom = acnv_data[0]
    acnv_start = acnv_data[1].split("-")[0]
    acnv_end = acnv_data[1].split("-")[1]

    overlapping_bedregions = []
    if acnv_chrom in bedfiledata:
        for bedregion in bedfiledata[acnv_chrom]:
            if acnv_start <= bedregion.exon_end and bedregion.exon_start <= acnv_end:
                overlapping_bedregions.append(bedregion)
    return overlapping_bedregions


def display_bedregions(bedregions, arraycnvregion):
    """Displays BED regions overlapping with the array CNV region.

    Parameters
    ----------
    bedregions : list of Exon
        BED regions overlapping with the array CNV region
    arraycnvregion : str
        Array CNV as a region (chr:start-end)
    """
    print(f"BED reegions overlapping with array CNV {arraycnvregion}")
    for bedregion in bedregions:
        print(f"{bedregion.exon_chrom}\t{bedregion.exon_start}\t{bedregion.exon_end}\t{bedregion.gene_name}")
