class ConradExon:
    def __init__(self, cexonchrom, cexonstart, cexonend, cexonname):
        self.exon_chrom = cexonchrom
        self.exon_start = cexonstart
        self.exon_end = cexonend
        self.exon_name = cexonname
        self.gene_name = cexonname.split("_")[0]
