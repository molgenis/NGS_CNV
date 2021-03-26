class ArrayCnv:
    def __init__(self, samplename, chrom, startpos, endpos, cnvcall, cnvsize, numprobes, numgenes, cnvclass):
        self.cnv_sample = samplename
        self.cnv_chrom = chrom
        self.cnv_start = startpos
        self.cnv_end = endpos
        self.cnv_call = cnvcall
        self.cnv_size = cnvsize
        self.number_of_probes = numprobes
        self.cnv_number_of_genes = numgenes
        self.cnv_class = cnvclass
        self.exons = []
        self.classification = ""
        self.call_result = ""
        self.wes_cnvs = []

    def get_length(self):
        return abs(self.cnv_end - self.cnv_start)

    def get_region(self):
        return f"{self.cnv_chrom}:{self.cnv_start}-{self.cnv_end}"

    def num_of_probes(self):
        return self.number_of_probes

    def num_of_exons(self):
        return len(self.exons)

    def get_gene_names(self):
        genelist = []
        for exon in self.exons:
            genelist.extend(exon.get_gene_names())
        return list(set(genelist))

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

    def num_of_wes_cnvs(self):
        return len(self.wes_cnvs)

    def __str__(self):
        return f"{self.cnv_sample}\t{self.cnv_chrom}\t{self.cnv_start}\t{self.cnv_end}"
