#!/usr/bin/env python
class UmcgCommonCnv:
    def __init__(self, chrom, startpos, endpos, nameid, score, strand, thickstart, thickend, rgbstr):
        self.cnvchrom = chrom
        self.cnvstart = startpos
        self.cnvend = endpos
        self.cnvname = nameid
        self.cnvscore = score
        self.strand = strand
        self.thickstart = thickstart
        self.thickend = thickend
        self.rgbvalue = rgbstr

    def get_region_str(self):
        return f"{self.chrom}:{self.startpos}-{self.endpos}"

    def is_male(self):
        return "_Male" in self.cnvname

    def is_female(self):
        return "_Female" in self.cnvname

    def __str__(self):
        return f"{self.cnvchrom}\t{self.cnvstart}\t{self.cnvend}\t{self.cnvname}\t{self.cnvscore}\t{self.strand}\t{self.thickstart}\t{self.thickend}\t{self.rgbvalue}"
