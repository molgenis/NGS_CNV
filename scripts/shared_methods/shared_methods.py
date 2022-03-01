#!/usr/bin/env python
def get_gatk_region(gatkregionstr):
    gatk_chrom = gatkregionstr.split(":")[0]
    gatk_start = int(gatkregionstr.split(":")[1].split("-")[0])
    gatk_end = int(gatkregionstr.split(":")[1].split("-")[1])
    return [gatk_chrom, gatk_start, gatk_end]
