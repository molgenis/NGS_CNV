class Probe:
    def __init__(self, chrom, startpos, endpos):
        self.probe_chrom = chrom
        self.probe_start = startpos
        self.probe_end = endpos

    def get_length(self):
        return abs(self.probe_end - self.probe_start)

    def __str__(self):
        return f"{self.probe_chrom}\t{self.probe_start}\t{self.probe_end}"
