class FpResult:
    def __init__(self, samplename, gatkcnv, gatkcall, gatksize, classification, exonsnum, probesnum):
        self.sample_name = samplename
        self.gatk_cnv = gatkcnv
        self.gatk_call = gatkcall
        self.gatk_size = gatksize
        self.classification = classification
        self.num_of_exons = exonsnum
        self.num_of_probes = probesnum

    def __str__(self):
        return f"{self.sample_name}\t{self.gatk_cnv}\t{self.gatk_call}\t{self.gatk_size}\t{self.classification}\t{self.num_of_exons}\t{self.num_of_probes}"
