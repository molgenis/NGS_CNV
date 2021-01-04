class Cnv:
    def __init__(self, samplename, samplepseudo, chrom, startpos, endpos, npcr, mlcr, call):
        self.cnv_sample = samplename
        self.cnv_sample_pseudo = samplepseudo
        self.cnv_chrom = f"chr{chrom}"
        self.cnv_start = startpos
        self.cnv_end = endpos
        self.num_points_copy_ratio = npcr
        self.median_log2_copy_ratio = mlcr
        self.cnv_call = call
        self.probes = []
        self.exons = []
        self.classification = []
        self.call_result = ""
        self.array_cnv = None
        self.left_hangover = None
        self.right_hangover = None

    def get_length(self):
        """Return the CNV length.
        
        Returns
        -------
        int
            Length of CNV
        """
        return abs(self.cnv_end - self.cnv_start)
    
    def get_region(self):
        return f"{self.cnv_chrom}:{self.cnv_start}-{self.cnv_end}"

    def num_of_probes(self):
        """Return the number of probes overlapping with the CNV.
        
        Returns
        -------
        int
            Number of probes
        """
        return len(self.probes)

    def num_of_exons(self):
        """Return the number of exons overlapping with the CNV.
        
        Returns
        -------
        int
            Number of exons
        """
        return len(self.exons)

    def get_gene_names(self):
        genelist = []
        for exon in self.exons:
            genelist.extend(exon.get_gene_names())
        return list(set(genelist))

    def segment_overlap(self, startpos, endpos):
        return self.cnv_start <= endpos and startpos <= self.cnv_end

    def cnv_overlap(self, other_cnv):
        """Determine and return whether the current CNV overlaps with a provided CNV.
        
        Returns
        -------
        bool
            True if CNVs overlap, False if not
        """
        return self.cnv_start <= other_cnv.cnv_end and other_cnv.cnv_start <= self.cnv_end

    def percentage_overlap(self, arraycnv_range):
        cnv_range = range(self.cnv_start, self.cnv_end)
        start_diff = arraycnv_range[0] - cnv_range[0]
        end_diff = arraycnv_range[-1] - cnv_range[-1]

        # Check whether the GATK CNV completely overlaps with the array CNV
        if start_diff > 0 and end_diff < 0:
                return 100

        # Check if the GATK CNV overlaps perfectly with the array CNV
        if start_diff == 0 and end_diff == 0:
                return 100

        # Determine the percentagr of overlap
        overlap_size = len(arraycnv_range)
        if start_diff < 0:
                overlap_size = overlap_size - abs(start_diff)
                if end_diff > 0:
                        overlap_size = overlap_size - end_diff

        overlap_perc = (overlap_size / len(arraycnv_range)) * 100
        return round(overlap_perc, 3)


    def __str__(self):
        """Return a string representation of the CNV entry.
        
        Returns
        -------
        str
            CNV as combined seg file entry
        """
        return f"{self.cnv_sample}\t{self.cnv_chrom}\t{self.cnv_start}\t{self.cnv_end}\t{self.cnv_num_points_copy_ratio}\t{self.cnv_median_log2_copy_ratio}\t{self.cnv_call}"
