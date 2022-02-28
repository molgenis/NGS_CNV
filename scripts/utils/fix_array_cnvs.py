#!/usr/bin/env python
def fix_array_cnv_notation(infileloc, outfileloc):
    try:
        infile = open(infileloc, 'r')
        outfile = open(outfileloc, 'w')
        
        outfile.write(next(infile))
        for inline in infile:
            inlinedata = inline.strip().split("\t")
            
            if inlinedata[4] != "NA":
                fixedregion = fix_array_region(inlinedata[4])
                linepart1 = "\t".join(inlinedata[0:4])
                linepart2 = "\t".join(inlinedata[5:])
                outfile.write(f"{linepart1}\t{fixedregion}\t{linepart2}\n")
            else:
                outfile.write(inline)
        infile.close()
        outfile.close()
    except IOError:
        print("Some error happened :(")


def fix_array_region(arrayregion):
    arraychrom = arrayregion.split("-")[0]
    arraystart = arrayregion.split("-")[1].split(":")[0]
    arrayend = arrayregion.split("-")[1].split(":")[1]
    return f"{arraychrom}:{arraystart}-{arrayend}"
