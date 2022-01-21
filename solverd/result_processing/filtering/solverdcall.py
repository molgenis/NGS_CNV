class SolveRdCall:
    def __init__(self, ename, cnvchrom, cnvstart, cnvend, cnvlen, cnvtype):
        self.sample_name = ename
        self.cnv_chrom = cnvchrom
        self.cnv_start = cnvstart
        self.cnv_end = cnvend
        self.cnv_length = cnvlen
        self.cnv_call = cnvtype
        self.overlapping_ccrs = []
        self.overlapping_ccrs_2 = []
        self.ern_genes = []

    def get_region_string(self):
        return f"{self.cnv_chrom}:{self.cnv_start}-{self.cnv_end}"

    def get_region_string_2(self):
        return self.cnv_chrom + ":" + "{:,}".format(self.cnv_start) + "-" + "{:,}".format(self.cnv_end)

    def get_ccrs_calls(self, usecommas=False):
        [x.get_region_string_2() for x in self.overlapping_ccrs] if usecommas else [x.get_region_string() for x in self.overlapping_ccrs]

    def get_ccrs_calls_2(self, usecommas=False):
        [x.get_region_string_2() for x in self.overlapping_ccrs_2] if usecommas else [x.get_region_string() for x in self.overlapping_ccrs_2]

# https://pythonguides.com/python-format-number-with-commas/
