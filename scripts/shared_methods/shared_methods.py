def get_gatk_region(gatkregionstr):
    gatk_chrom = gatkregionstr.split(":")[0]
    gatk_start = int(gatkregionstr.split(":")[1].strip("-")[0])
    gatk_end = int(gatkregionstr.split(":")[1].strip("-")[1])
    return [gatk_chrom, gatk_start, gatk_end]
