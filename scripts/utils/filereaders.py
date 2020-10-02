def read_interval_data(intervalfileloc):
    """Read and return interval data from a provided file.
    
    Interval data should be data in a denoised of standardized tsv interval file.
    
    Parameters
    ----------
    intervalfileloc : str
        Path to interval file.
    
    Returns
    -------
    intervaldata : dict
        Read interval data.
    """
    intervaldata = {}
    try:
        with open(intervalfileloc, "r") as intervalfile:
            for fileline in intervalfile:
                if not fileline.startswith(("@", "CONTIG")):
                    filelinedata = fileline.strip().split("\t")
                    intervalkey = f"{filelinedata[0]}_{filelinedata[1]}"
                    intervaldata[intervalkey] = filelinedata
    except IOError:
        print("Could not read interval file.")
    finally:
        return intervaldata


def read_allelic_data(allelicfileloc):
    """Read and return allelic data from a provided file.
    
    Parameters
    ----------
    allelicfileloc : str
        Path to allelic data file.
    
    Returns
    -------
    allelicdata : dict
        Read allelic data.
    """
    allelicdata = {}
    lineindex = 1
    try:
        with open(allelicfileloc, "r") as allelicfile:
            for fileline in allelicfile:
                if not fileline.startswith(("@", "CONTIG")):
                    filelinedata = fileline.strip().split("\t")
                    allelicdata[lineindex] = filelinedata
                    lineindex += 1
    except IOError:
        print("Could not read allelic file.")
    finally:
        return allelicdata


def read_classification_file(infileloc):
    gatkresultdata = {}
    try:
        with open(infileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")

                if filelinedata[0] not in gatkresultdata:
                    gatkresultdata[filelinedata[0]] = []

                gatkreg = get_gatk_region(filelinedata[1])
                gatkgenes = filelinedata[13].split(":")
                gatkcall = GatkCall(filelinedata[0], gatkreg[0], gatkreg[1], gatkreg[2], filelinedata[2], int(filelinedata[3]), filelinedata[4], filelinedata[9], filelinedata[10], int(filelinedata[11]), int(filelinedata[12]), gatkgenes, fileline.strip())
                gatkresultdata[filelinedata[0]].append(gatkcall)
    except IOError:
        print(f"Could not read input file {infileloc}")
    finally:
        return gatkresultdata


def read_array_cnvs(arrayfileloc):
    array_cnv_data = {}
    try:
        with open(arrayfileloc, 'r') as arrayfile:
            next(arrayfile)
            for arrayline in arrayfile:
                arraylinedata = arrayline.strip().split("\t")

                if arraylinedata[1] != "nan":
                    if arraylinedata[0] not in array_cnv_data:
                        array_cnv_data[arraylinedata[0]] = {}

                    arrayregion = get_gatk_region(arraylinedata[1])
                    regionstr = f"{arrayregion[0]}:{arrayregion[1]}-{arrayregion[2]}"

                    if regionstr not in array_cnv_data[arraylinedata[0]]:
                        arraycnv = ArrayCnv(arraylinedata[0], arrayregion[0], arrayregion[1], arrayregion[2], arraylinedata[2], int(arraylinedata[3]), int(arraylinedata[4]), int(arraylinedata[5]), arraylinedata[6])
                        array_cnv_data[arraylinedata[0]][regionstr] = arraycnv
    except IOError:
        print(f"Could not open {arrayfileloc}")
    finally:
        return array_cnv_data


def read_umcg_common_cnv_file(ccnvfileloc):
    ccnv_data = {}
    try:
        with open(ccnvfileloc, 'r') as ccnvfile:
            next(ccnvfile)
            for ccnvline in ccnvfile:
                ccnvlinedata = ccnvline.strip().split("\t")
                if ccnvlinedata[0] not in ccnv_data:
                    ccnv_data[ccnvlinedata[0]] = []
                ccnv_data[ccnvlinedata[0]].append(UmcgCommonCnv(ccnvlinedata[0], int(ccnvlinedata[1]), int(ccnvlinedata[2]), ccnvlinedata[3], int(ccnvlinedata[4]), ccnvlinedata[5], int(ccnvlinedata[6]), int(ccnvlinedata[7]), ccnvlinedata[8]))
    except IOError:
        print(f"Could not read common cnv file: {ccnvfileloc}")
    finally:
        return ccnv_data
