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

    def num_of_genes(self):
        return len(self.gene_names)

    def get_region_str(self):
        return f"{self.chrom}:{self.startpos}-{self.endpos}"
