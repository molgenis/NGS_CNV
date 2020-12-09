class ExomeDepthCall:
    def __init__(self, samplename, pseudosample, cnvstartp, cnvendp, cnvcall, cnvexonnum, cnvstart, cnvend, cnvchrom, cnvid, cnvbf, cnvreadexp, cnvreadobs, cnvreadratio, cnvconrad):
        self.cnv_sample = samplename
        self.cnv_sample_pseudo = pseudosample
        self.cnv_chrom = cnvchrom
        self.cnv_start = cnvstart
        self.cnv_end = cnvend
        self.cnv_call = cnvcall
        self.probes = []
        self.exons = []
        self.classification = ""
        self.call_result = ""
        self.array_cnv = None
        self.left_hangover = None
        self.right_hangover = None
        self.percent_overlap = None
        self.startp = cnvstartp
        self.endp = cnvendp
        self.num_exons = cnvexonnum
        self.identifier = cnvid
        self.bf = cnvbf
        self.reads_expected = cnvreadenv
        self.reads_observed = cnvreadobs
        self.reads_ratio = cnvreadratio
        self.conrad_hg19 = cnvconrad


    def get_region(self):
        return f"{self.chrom}:{self.startpos}-{self.endpos}"

    def __str__(self):
        return f"{self.startp}\t{self.endp}\t{self.call_type}\t{self.num_exons}\t{self.startpos}\t{self.endpos}\t{self.chrom}\t{self.identifier}\t"
               f"{self.bf}\t{self.reads_expected}\t{self.reads_observed}\t{self.reads_ratio}\t{self.conrad_hg19}"
