class FpSummary:
    def __init__(self, classification, callsnum, lenavg, lenmed, exonavg, exonmed, probeavg, probemed):
        self.classification = classification
        self.num_of_calls = callsnum
        self.avg_len = lenavg
        self.med_len = lenmed
        self.avg_exons = exonavg
        self.med_exons = exonmed
        self.avg_probes = probeavg
        self.med_probes = probemed

    def __str__(self):
        return f"{self.classification}\t{self.num_of_calls}\t{self.avg_len}\t{self.med_len}\t{self.avg_exons}\t{self.med_exons}\t{self.avg_probes}\t{self.med_probes}"
