#!/usr/bin/env python
import sys


master_path = "/groups/umcg-solve-rd/rsc01/releases/master/bam/"

try:
    infile = open(sys.argv[1], 'r')
    outfile = open(sys.argv[2], 'w')
    for fileline in infile:
        if fileline.strip().endswith((".bam", ".bam.bai", ".bai")):
            bamfile = fileline.strip().split()[-1]
            outfile.write(f"{master_path}{bamfile}\n")
    infile.close()
    outfile.close()
except IOError:
    print("AAP")
