class ConiferCall:
    def __init__(self, samplename, cnvchrom, cnvstart, cnvend, callresult):
        self.sampleid = samplename
        self.chrom = cnvchrom
        self.startpos = cnvstart
        self.endpos = cnvend
        self.cnvcall = callresult

    def get_region(self):
        """Returns the Conifer call CNV as a region string.

        Returns
        -------
        str
            CNV region as chr:start-end
        """
        return f"{self.chrom}:{self.startpos}-{self.endpos}"

    def __str__(self):
        return f"{self.sampleid}\t{self.chrom}\t{self.startpos}\t{self.endpos}"
