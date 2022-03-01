#!/usr/bin/env python
class CcrsEntry:
    def __init__(self, cnv_chrom, cnv_start, cnv_end, cnv_npcr, cnv_mlcr):
        self.chrom = cnv_chrom
        self.startpos = cnv_start
        self.endpos = cnv_end
        self.num_points_copy_ratio = cnv_npcr
        self.mean_log2_copy_ratio = cnv_mlcr

    def get_size(self):
        return self.endpos - self.startpos

    def get_region(self):
        return f"{self.chrom}:{self.startpos}-{self.endpos}"

    def __str__(self):
        return f"{self.chrom}\t{self.startpos}\t{self.endpos}\t{self.num_points_copy_ratio}\t{self.mean_log2_copy_ratio}"
