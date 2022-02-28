#!/usr/bin/env python
class Exon:
    def __init__(self, chrom, startpos, endpos, genename):
        self.exon_chrom = chrom
        self.exon_start = startpos
        self.exon_end = endpos
        self.gene_name = genename

    def get_length(self):
        return abs(self.exon_end - self.exon_start)

    def get_gene_names(self):
        return self.gene_name.split(":")

    def __str__(self):
        return f"{self.exon_chrom}\t{self.exon_start}\t{self.exon_end}\t{self.gene_name}"
