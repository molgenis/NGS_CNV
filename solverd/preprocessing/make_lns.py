import sys

linkdir = sys.argv[3]
infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[2], 'w')
for fileline in infile:
    bamfile = fileline.strip().split("/")[-1]
    outfile.write(f"ln -s {fileline.strip()} {linkdir}{bamfile}\n")
infile.close()
outfile.close()
