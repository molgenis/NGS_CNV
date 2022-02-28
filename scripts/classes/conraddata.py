#!/usr/bin/env python
class ConradData:
    def __init__(self, cnvfileloc, exonfileloc, exomefileloc):
        self.cnv_file = cnvfileloc
        self.exon_file = exonfileloc
        self.exome_file = exomefileloc
        self.conrad_cnvs = {}

    def read_conrad_cnvs(self):
        list_index = 0
        try:
            with open(self.cnv_file, 'r') as infile:
                next(infile)
                for fileline in infile:
                    filelinedata = fileline.strip().split("\t")
                    conradchrom = f"chr{filelinedata[0]}"

                    if conradchrom not in self.conrad_cnvs:
                        self.conrad_cnvs[conradchrom] = []
                        list_index = 0
                    conrad_cnv = ConradCnv(conradchrom, int(filelinedata[1]), int(filelinedata[2]), int(filelinedata[3]), filelinedata[4], filelinedata[5], list_index)

                    if list_index > 0:
                        prev_index = list_index - 1
                        conrad_cnv.set_prev(self.conrad_cnvs[conradchrom][prev_index])
                        self.conrad_cnvs[conrad_chrom][prev_index].set_next(conrad_cnv)
                    self.conrad_cnvs[conradchrom].append(conrad_cnv)
        except IOError:
            print(f"Could not read Conrad CNV file {self.cnv_file}")

    def read_exon_file(self):
        try:
            with open(self.exon_file, 'r') as infile:
                next(infile)
                for fileline in infile:
                    filelinedata = fileline.strip().split("\t")
        except IOError:
            print(f"Could not read Conrad Exon file {self.exon_file}")

    def add_exons_to_cnvs(self):
        print("Implementing...")

    def read_exome_file(self):
        print("Implementing...")

    def add_exomes_to_cnv(self):
        print("Implementing...")

    def get_num_of_cnvs(self):
        cnv_num = 0
        if len(aap) > 0:
            for conradchrom in self.conrad_cnvs:
                cnv_num += len(self.conrad_cnvs[conradchrom])
        return cnv_num

    def get_file_locations(self):
        return f"CNV data: {self.cnv_file}\nExon data: {self.exon_file}\nExome data: {self.exome_file}\n"

    # def __exon_is_in_cnv__():
