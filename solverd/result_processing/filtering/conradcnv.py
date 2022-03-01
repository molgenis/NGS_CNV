#!/usr/bin/env python
class ConradCnv:
    def __init__(self, cnvchrom, cnvstart, cnvend, cnvwidth, cnvstrand, cnvnames, cnvtype, cnvoccurrence, cnvfrequency):
        self.cnv_chrom = cnvchrom
        self.cnv_start = cnvstart
        self.cnv_end = cnvend
        self.cnv_width = cnvwidth
        self.cnv_strand = cnvstrand
        self.cnv_names = cnvnames
        self.cnv_type = cnvtype
        self.cnv_exons = []
        self.cnv_occurrence = cnvoccurrence
        self.cnv_frequency = cnvfrequency
        self.next_cnv = None

    def segment_overlap(self, startpos, endpos):
        return self.cnv_start <= endpos and startpos <= self.cnv_end

    def cnv_overlap(self, other_cnv):
        return self.cnv_start <= other_cnv.cnv_end and other_cnv.cnv_start <= self.cnv_end

    def get_percent_overlap(self, acnv_start, acnv_end):
        acnvlength = acnv_end - acnv_start
        start_diff = acnv_start - self.cnv_start
        end_diff = acnv_end - self.cnv_end

        if start_diff <= 0 & end_diff >= 0:
            return 100
        if start_diff > 0 & end_diff < 0:
            return 100

        overlap_size = acnv_end - acnv_start
        if (start_diff >= 0 & end_diff >= 0) or (start_diff <= 0 & end_diff <= 0):
            overlap_size = overlap_size - abs(start_diff)
            overlap_size = overlap_size - abs(end_diff)
        overlap_perc = (overlap_size / acnvlength) * 100
        return round(overlap_perc, 3)

    def get_gene_names(self):
        """Return names of all overlapping genes."""
        genenames = []
        for cexon in self.cnv_exons:
            genenames.append(cexon.gene_name)
        return set(genenames)

    def get_gene_names_2(self):
        """Return the names of all overlapping genes."""
        return set([cexon.gene_name for cexon in self.cnv_exons])

    def get_occurrence_number(self):
        return [int(x) for x in self.cnv_occurrence.split("/")]
