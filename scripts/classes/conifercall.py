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
        return f"{self.cnv_chrom}:{self.cnv_start}-{self.cnv_end}"

    def get_length(self):
        return abs(self.cnv_end - self.cnv_start)

    def segment_overlap(self, startpos, endpos):
        return self.cnv_start <= endpos and startpos <= self.cnv_end

    def cnv_overlap(self, othercnv):
        if self.cnv_chrom == othercnv.cnv_chrom:
            return self.cnv_start <= othercnv.cnv_end and othercnv.cnv_start <= self.cnv_end
        return False

    def get_percent_overlap(self, arraycnv_range):
        cnv_range = range(self.cnv_start, self.cnv_end)
        start_diff = arraycnv_range[0] - cnv_range[0]
        end_diff = arraycnv_range[-1] - cnv_range[-1]

        if start_diff > 0 and end_diff < 0:
                return 100
        if start_diff == 0 and end_diff == 0:
                return 100

        overlap_size = len(arraycnv_range)
        if start_diff < 0:
                overlap_size = overlap_size - abs(start_diff)
                if end_diff > 0:
                        overlap_size = overlap_size - end_diff
        overlap_perc = (overlap_size / len(arraycnv_range)) * 100
        return round(overlap_perc, 3)

    def get_gene_names(self):
        genelist = []
        for exon in self.exons:
            genelist.extend(exon.get_gene_names())
        return list(set(genelist))

    def __str__(self):
        return f"{self.cnv_sample}\t{self.cnv_chrom}\t{self.cnv_start}\t{self.cnv_end}"
