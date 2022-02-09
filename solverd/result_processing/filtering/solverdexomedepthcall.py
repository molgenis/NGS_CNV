class SolveRdExomeDepthCall:
    def __init__(self, ename, cnvchrom, cnvstart, cnvend, cnvlen, cnvtype, cnvnexons):
        self.sample_name = ename
        self.cnv_chrom = cnvchrom
        self.cnv_start = cnvstart
        self.cnv_end = cnvend
        self.cnv_length = cnvlen
        self.cnv_call = cnvtype
        self.cnv_nexons = cnvnexons
        self.overlapping_ccrs = []
        self.overlapping_ccrs_2 = []
        self.ern_genes = []

    def get_region_string(self):
        """Return the Solve-RD call as a genomic region."""
        return f"{self.cnv_chrom}:{self.cnv_start}-{self.cnv_end}"

    def get_region_string_2(self):
        """Return the Solve-RD calls as a genomic region with comma formatting."""
        return self.cnv_chrom + ":" + "{:,}".format(self.cnv_start) + "-" + "{:,}".format(self.cnv_end)

    def get_ccrs_calls(self, usecommas=False):
        """Return all proper (same call dup/del) overlapping CCRS calls as comma formatted region strings."""
        [x.get_region_string_2() for x in self.overlapping_ccrs] if usecommas else [x.get_region_string() for x in self.overlapping_ccrs]

    def get_ccrs_calls_2(self, usecommas=False):
        """Return all improper (different call) overlapping CCRS calls aas comma formatted region strings."""
        [x.get_region_string_2() for x in self.overlapping_ccrs_2] if usecommas else [x.get_region_string() for x in self.overlapping_ccrs_2]
