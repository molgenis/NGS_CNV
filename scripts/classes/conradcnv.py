#!/usr/bin/env python
class ConradCnv:
    def __init__(self, cchrom, cstart, cend, cwidth, cstrand, cname, lindex):
        self.conrad_chrom = cchrom
        self.conrad_start = cstart
        self.conrad_end = cend
        self.conrad_width = cwidth
        self.conrad_strand = cstrand
        self.conrad_name = cname
        self.list_index = lindex
        self.prev_cnv = None
        self.next_cnv = None
        self.gatk_cnv = None
        self.conrad_exons = []

    def set_prev(self, prevcnv):
        self.prev_cnv = prevcnv

    def set_next(self, nextcnv):
        self.next_cnv = nextcnv

    def set_gatk_cnv(self, gatkcnv):
        self.gatk_cnv = gatkcnv

    def has_prev(self):
        return self.prev_cnv is not None

    def has_next(self):
        return self.next_cnv is not None

    def has_gatk_cnv(self):
        return self.gatk_cnv is not None

    def get_prev(self):
        return self.prev_cnv

    def get_next(self):
        return self.next_cnv

    def get_gatk_cnv(self):
        return self.gatk_cnv

    def add_conrad_exon(self, conradexon):
        self.conrad_exons.append(conradexon)

    def num_of_exons(self):
        return len(self.conrad_exons)

    def get_genenames(self):
        genenames = []
        for cexon in self.conrad_exons:
            genenames.extend(cexon.get_gene_names())
        return list(set(genenames))

    def __str__(self):
        return f"{self.conrad_chrom}\t{self.conrad_start}\t{self.conrad_end}\t{self.conrad_width}\t{self.conrad_strand}\t{self.conrad_name}"
