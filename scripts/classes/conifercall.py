class ConiferCall:
    def __init__(self, samplepseudo, samplename, cnvchrom, cnvstart, cnvend, callresult):
        self.cnv_sample = samplename
        self.cnv_sample_pseudo = samplepseudo
        self.cnv_chrom = cnvchrom
        self.cnv_start = cnvstart
        self.cnv_end = cnvend
        self.cnv_call = callresult
        self.probes = []
        self.exons = []
        self.classification = ""
        self.call_result = ""
        self.array_cnv = None
        self.left_hangover = None
        self.right_hangover = None
        self.percent_overlap = None

    def get_region(self):
        """Returns the Conifer call CNV as a region string.

        Returns
        -------
        str
            CNV region as chr:start-end
        """
        return f"{self.chrom}:{self.startpos}-{self.endpos}"

    def segment_overlap(self, startpos, endpos):
        return self.cnv_start <= endpos and startpos <= self.cnv_end

    def cnv_overlap(self, cnvchrom, cnvstart, cnvend):
        if self.cnv_chrom == cnvchrom:
            return self.cnv_start <= cnvend and cnvstart <= self.cnv_end
        return False

    def get_gene_names(self):
        genelist = []
        for exon in self.exons:
            genelist.extend(exon.get_gene_names())
        return list(set(genelist))

    def __str__(self):
        return f"{self.cnv_sample}\t{self.cnv_chrom}\t{self.cnv_start}\t{self.cnv_end}"
