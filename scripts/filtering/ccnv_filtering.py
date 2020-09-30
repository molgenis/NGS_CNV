def filter_classification_results(infileloc, ccnv_data, outfileloc):
    try:
        outfile = open(outfileloc, 'w')
        with open(infileloc, 'r') as infile:
            outfile.write(next(infile))
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                gatkregion = get_gatk_region(filelinedata[1])
                
                if gatkregion[0] in ccnv_data:
                    overlapcount = count_ccnv_overlaps(gatkregion, ccnv_data[gatkregion[0]])
                    if overlapcount < 10:
                        outfile.write(fileline)
        outfile.close()
    except IOError:
        print("Something went wrong :(")


def get_gatk_region(gatkregionstr):
    gatk_chrom = gatkregionstr.split(":")[0]
    gatk_start = int(gatkregionstr.split(":")[1].strip("-")[0])
    gatk_end = int(gatkregionstr.split(":")[1].strip("-")[1])
    return [gatk_chrom, gatk_start, gatk_end]


def count_ccnv_overlaps(gatkregion, ccnvs):
    overlap_num = 0
    for ccnv in ccnvs:
        if ccnv_overlaps(gatkregion[1], gatkregion[2], ccnv.cnvstart, ccnv.cnvend):
            overlap_num += 1
    return overlap_num


def ccnv_overlaps(gatkstart, gatkend, ccnvstart, ccnvend):
    return gatkstart <= ccnvend and ccnvstart <= gatkend
