import statistics

class CcrsCall:
    def __init__(self, ccrssample, ccrschrom, ccrsstart, ccrsend, ccrsprobes, ccrscall, ccrssegmean):
        self.ccrs_sample = ccrssample
        self.ccrs_chrom = ccrschrom
        self.ccrs_start = ccrsstart
        self.ccrs_end = ccrsend
        self.ccrs_numofprobes = ccrsprobes
        self.ccrs_call = ccrscall
        self.ccrs_segmentmean = ccrssegmean
        self.ccrs_occurrence = ""
        self.ccrs_frequency = -1.0
        self.ccrs_callgroup_name = ""
        self.ccrs_callgroup_occurrence = ""
        self.ccrs_callgroup_frequency = -1.0
        self.conrad_occurrence = []
        self.conrad_frequency = []
        self.gnomad_frequency = -1.0
        self.conrad_cnvs = []
        self.conrad_exons = []
        self.keep_call = False
        self.ern_genes = []
        self.gnomad_entries = []
        self.callgroup_processed = False
        self.callgroup_common = False

    def get_call_length(self):
        return self.ccrs_end - self.ccrs_start

    def get_call_translation(self):
        """Return word translation for the type of call."""
        if self.ccrs_call == "+":
            return "Duplication"
        if self.ccrs_call == "-":
            return "Deletion"
        if self.ccrs_call == "0":
            return "Neutral"

    def get_gene_names(self):
        """Return names of all overlapping Conrad genes."""
        genenames = []
        for cexon in self.conrad_exons:
            genenames.append(cexon.gene_name)
        return set(genenames)

    def get_gene_names_2(self):
        """Return the names of all overlapping Conrad genes."""
        return set([cexon.gene_name for cexon in self.conrad_exons])

    def get_conrad_gene_names(self):
        """Return gene names of all overlapping Conrad CNVs."""
        conrad_genenames = []
        for conradcnv in self.conrad_cnvs:
            conrad_genenames.extend(conradcnv.get_gene_names_2())
        return conrad_genenames

    def get_region_string(self):
        """Return the CCRS call as a genomic region string."""
        return f"{self.ccrs_chrom}:{self.ccrs_start}-{self.ccrs_end}"

    def get_region_string_2(self):
        return self.ccrs_chrom + ":" + "{:,}".format(self.ccrs_start) + "-" + "{:,}".format(self.ccrs_end)

    def get_conrad_cnv_strings(self):
        """Return all overlapping Conrad CNVs as string representations"""
        if len(self.conrad_cnvs) > 0:
            return [f"{x.cnv_chrom}:{x.cnv_start}-{x.cnv_end}" for x in self.conrad_cnvs]
        return []

    def get_conrad_call_types(self):
        """Return the call types of all overlapping Conrad CNVs."""
        calltypes = []
        for ccnv in self.cnorad_cnvs:
            calltypes.extend(ccnv.cnv_type)
        return calltypes

    def can_be_conrad_filtered(self, calltranslationtable):
        """Determine whether the CCRS call can be filtered by Conrad (call types, gene names)."""
        return self.conradcnvs_have_same_calls(calltranslationtable) and self.all_genes_in_conradcnvs()

    def conradcnvs_have_same_calls(self, calltranslationtable):
        """Return whether the overlapping Conrad CNVs have the same call as the CCRS call."""
        for ccnv in self.conrad_cnvs:
            if calltranslationtable[self.ccrs_call] not in ccnv.cnv_type:
                return False
        return True

    def all_genes_in_conradcnvs(self):
        """Return whether all CCRS call gene names are in the overlapping Conrad CNVs."""
        return len(set(self.get_gene_names_2()) - set(self.get_conrad_gene_names())) == 0

    def get_ern_genes(self):
        """Return a list with names of all overlapping ERN genes."""
        return [x.ern_name for x in self.ern_genes]

    def to_ccrs_file_line(self):
        """Return the CCRS call as a combined CCRS file line representation."""
        ccrs_file_line = f"{self.ccrs_sample}\t{self.ccrs_chrom}\t{self.ccrs_start}\t{self.ccrs_end}\t{self.ccrs_numofprobes}\t{self.ccrs_call}\t{self.ccrs_segmentmean}"
        if self.ccrs_occurrence != "" and self.ccrs_frequency != -1.0:
            ccrs_file_line += f"\t{self.ccrs_occurrence}\t{self.ccrs_frequency}"
        if len(self.conrad_occurrence) > 0 and len(self.conrad_frequency) > 0:
            ccrs_file_line += "\t"
            ccrs_file_line += "|".join(self.conrad_occurrence)
            ccrs_file_line += "\t"
            ccrs_file_line += "|".join(str(x) for x in self.conrad_frequency)
        if self.gnomad_frequency != -1.0:
            ccrs_file_line += f"\t{self.gnomad_frequency}"
        ccrs_file_line += "\n"
        return ccrs_file_line

    def to_ccrs_file_line_2(self):
        """Return the CCRS call as a combined CCRS file line representation."""
        ccrs_file_line = f"{self.ccrs_sample}\t{self.ccrs_chrom}\t{self.ccrs_start}\t{self.ccrs_end}\t{self.ccrs_numofprobes}\t{self.ccrs_call}\t{self.ccrs_segmentmean}"
        if self.ccrs_occurrence != "" and self.ccrs_frequency != -1.0:
            ccrs_file_line += f"\t{self.ccrs_occurrence}\t{self.ccrs_frequency}"
        if self.ccrs_callgroup_name != "" and self.ccrs_callgroup_occurrence != "" and self.ccrs_callgroup_frequency != -1.0:
            ccrs_file_line += f"\t{self.ccrs_callgroup_name}\t{self.ccrs_callgroup_occurrence}\t{self.ccrs_callgroup_frequency}"
        if len(self.conrad_occurrence) > 0 and len(self.conrad_frequency) > 0:
            ccrs_file_line += "\t"
            ccrs_file_line += "|".join(self.conrad_occurrence)
            ccrs_file_line += "\t"
            ccrs_file_line += "|".join(str(x) for x in self.conrad_frequency)
        if self.gnomad_frequency != -1.0:
            ccrs_file_line += f"\t{self.gnomad_frequency}"
        ccrs_file_line += "\n"
        return ccrs_file_line

    def get_conrad_occurrences(self):
        """Return the occurrences from the overlapping ConradCnv (self.conrad_cnvs)."""
        if len(self.conrad_cnvs) > 0:
            if len(self.conrad_cnvs) > 1:
                return "|".join([x.cnv_occurrence for x in self.conrad_cnvs])
            else:
                return self.conrad_cnvs[0].cnv_occurrence
        return "NA"

    def get_mean_conrad_occurrence(self):
        """Return the mean Conrad occurrence."""
        if len(self.conrad_occurrence) > 0:
            if len(self.conrad_occurrence) > 1:
                return statistics.mean(self.conrad_occurrence)
            else:
                return self.conrad_occurrence[0]
        return "NA"

    def determine_mean_conrad_occurrence(self):
        """Determine and return the mean occurrence of the overlapping Conrad CNVs."""
        if len(self.conrad_cnvs) > 0:
            if len(self.conrad_cnvs) > 1:
                nsamples = self.conrad_cnvs[0].cnv_occurrence.split("/")[1]
                mean_occurrence = statistics.mean([int(x.cnv_occurrence.split("/")[0]) for x in self.conrad_cnvs])
                return f"{mean_occurrence}/{nsamples}"
            else:
                return self.conrad_cnvs[0].cnv_occurrence
        return "NA"

    def get_conrad_frequencies(self):
        """Return the frequencies from the overlapping ConradCnv (self.conrad_cnvs)."""
        if len(self.conrad_cnvs) > 0:
            if len(self.conrad_cnvs) > 1:
                return "|".join([str(x.cnv_frequency) for x in self.conrad_cnvs])
            else:
                return str(self.conrad_cnvs[0].cnv_frequency)
        return "NA"

    def get_mean_conrad_frequency(self):
        """Return the mean Conrad frequency."""
        if len(self.conrad_frequency) > 0:
            if len(self.conrad_frequency) > 1:
                return statistics.mean(self.conrad_frequency)
            else:
                return self.conrad_frequency[0]
        return "NA"

    def determine_mean_conrad_frequency(self):
        """Determine and return the mean frequency overlapping Conrad CNVs."""
        if len(self.conrad_cnvs) > 0:
            if len(self.conrad_cnvs) > 1:
                return statistics.mean([x.cnv_frequency for x in self.conrad_cnvs])
            else:
                return self.conrad_cnvs[0].cnv_frequency
        return "NA"

    def get_summed_conrad_occurrence_and_frequency(self):
        """Return the summed conrad_frequency."""
        total_count = sum([int(x.cnv_occurrence.split("/")[0]) for x in self.conrad_cnvs])
        total_samples = int(self.conrad_cnvs[0].cnv_occurrence.split("/")[1])
        return [f"{total_count}/{total_samples}", round((total_count/total_samples)*100, 2)]

    def get_max_conrad_occurrence_and_frequency(self):
        """Return the maximum conrad frequency."""
        max_count = max([int(x.cnv_occurrence.split("/")[0]) for x in self.conrad_cnvs])
        total_samples = int(self.conrad_cnvs[0].cnv_occurrence.split("/")[1])
        return [f"{max_count}/{total_samples}", round((max_count/total_samples)*100, 2)]

    def get_conrad_occurrence_and_frequency(self, calltranslationstabel):
        """Return the Conrad occurrence and frequency."""
        if len(self.conrad_cnvs) > 0:
            c_occur = 0
            for ccnv in self.conrad_cnvs:
                if calltranslationstabel[self.ccrs_call] in ccnv.cnv_type:
                    if len(set(self.get_gene_names_2()) - set(ccnv.get_gene_names_2())) == 0:
                        c_occur += ccnv.get_occurrence_number()[0]
            return [f"{c_occur}/40", round((c_occur/40)*100, 2)]
        return ["0/40", 0.00]
