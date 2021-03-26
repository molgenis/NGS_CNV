class GatkCall:
    def __init__(self, samplename, gchrom, gstart, gend, gcall, gsize, acnv, gcallres, gclass, nexon, nprobe, genenames, callline):
        self.samplename = samplename
        self.chrom = gchrom
        self.startpos = gstart
        self.endpos = gend
        self.cnvcall = gcall
        self.cnvsize = gsize
        self.arraycnv = acnv
        self.callresult = gcallres
        self.classification = gclass
        self.exonnum = nexon
        self.probenum = nprobe
        self.gene_names = genenames
        self.conrad_cnvs = []
        self.call_line = callline

    def add_conrad_cnv(self, conradcnv):
        self.conrad_cnvs.append(conradcnv)

    def num_of_conrad_cnvs(self):
        return len(self.conrad_cnvs)

    def get_conrad_cnv_genes(self):
        conrad_genes = []
        for ccnv in self.conrad_cnvs:
            conrad_genes.extend(ccnv.get_genenames())
        return list(set(conrad_genes))

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


    def num_of_genes(self):
        return len(self.gene_names)

    def get_region_str(self):
        return f"{self.chrom}:{self.startpos}-{self.endpos}"
