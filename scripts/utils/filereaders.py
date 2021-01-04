from classes.arraycnv import ArrayCnv
from classes.exon import Exon
from classes.probe import Probe


def read_sample_table(samplefileloc):
    """Read and return a sample translation table.

    Parameters
    ----------
    samplefileloc : str
        Path to sample translation table file
    """
    sample_table = {}
    try:
        with open(samplefileloc, 'r') as samplefile:
            next(samplefile)
            for sampleline in samplefile:
                sampledata = sampleline.strip().split("\t")
                full_sample_id = sampledata[2].split(".")[0]

                if full_sample_id not in sample_table:
                    sample_table[full_sample_id] = sampledata[0]
    except IOError:
        print(f"Could not open sample table {samplefileloc}")
    finally:
        return sample_table


def read_probes_data(probefileloc):
    """Read and return probe data.

    The probe file should have three columns: chrom, start, end.

    Parameters
    ----------
    probefileloc : str
        Path to probes file
    """
    probe_data = {}
    try:
        with open(probefileloc, 'r') as probefile:
            for probeline in probefile:
                probelinedata = probeline.strip().split("\t")

                if probelinedata[0] not in probe_data:
                    probe_data[probelinedata[0]] = []
                probe_data[probelinedata[0]].append(Probe(probelinedata[0], int(probelinedata[1]), int(probelinedata[2])))
    except IOError:
        print(f"Could not open probe file {probefileloc}")
    finally:
        return probe_data


def read_exon_data(exonfileloc):
    """Read and return exon data from a file.

    The exon file should have four columns: chr, start, end, annotation/gene.

    Parameters
    ----------
    exonfileloc : str
        Path to exon file

    Returns
    -------
    exon_data : dict
        Exon data per chromosome
    """
    exon_data = {}
    try:
        with open(exonfileloc, 'r') as exonfile:
            for exonline in exonfile:
                exonlinedata = exonline.strip().split()
                exonchrom = exonlinedata[0]
                if not exonchrom.startswith("chr"):
                    exonchrom = f"chr{exonchrom}"

                if exonchrom not in exon_data:
                    exon_data[exonchrom] = []
                exon_data[exonchrom].append(Exon(exonlinedata[0], int(exonlinedata[1]), int(exonlinedata[2]), exonlinedata[3]))
    except IOError:
        print(f"Could not open exon file {exonfileloc}")
    finally:
        return exon_data


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
    """Read file containing Common CNVs.

    Parameters
    ----------
    ccnvfileloc : str

    Returns
    -------
    ccnv_data : dict
        Common CNV data
    """
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


def read_fp_classifications(classifications_fileloc):
    fp_data = []
    try:
        with open(classifications_fileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[10] == "FALSE POSITIVE":
                    fp_data.append(filelinedata[1])
    except IOError:
        print("Could not read supplied classifications file")
    finally:
        return fp_data


def read_fp_classifications_2(classfications_fileloc):
    fp_data = {}
    try:
        with open(classfications_fileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[10] == "FALSE POSITIVE":
                    if filelinedata[1] not in fp_data:
                        fp_data[filelinedata[1]] = []
                    fp_data[filelinedata[1]].append(filelinedata[0])
    except IOError:
        print("Could not read supplied classifications file")
    finally:
        return fp_data


def read_conifer_calls(coniferfileloc):
    """Read a file containing Conifer CNV calls.

    Parameters
    ----------
    coniferfileloc : str
        Path to conifer CNV calls file

    Returns
    -------
    conifer_calls : dict
        Conifer CNV calls per sample
    """
    conifer_calls = {}
    try:
        with open(coniferfileloc, 'r') as coniferfile:
            next(coniferfile)
            for fileline in coniferfile:
                filelinedata = fileline.strip().split("\t")

                if filelinedata[0] not in conifer_calls:
                    confcall = ConiferCall(filelinedata[0], filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), filelinedata[4])
                    conifer_calls[filelinedata[0]] = []
                conifer_calls[filelinedata[0]].append(confcall)
    except IOError:
        print(f"Could not read conifer file {coniferfileloc}")
    finally:
        return conifer_calls


def read_exomedepth_calls(exomedepthfileloc):
    """Read an ExomeDepth calls file of a single sample.

    Parameters
    ----------
    exomedepthfileloc : str
        Path to ExomeDepth sample calls.

    Returns
    -------
    exomedepth_calls : list of ExomeDepthCall
        ExomeDepth CNV calls for a single sample
    """
    exomedepth_calls = []
    try:
        with open(exomedepthfileloc, 'r') as exomedepthfile:
            next(exomedepthfile)

            for fileline in exomedepthfile:
                fileline = fileline.strip().replace("\"", "")
                filelinedata = fileline.split(",")

                conrad_data = filelinedata[12:]
                edc = ExomeDepthCall(int(filelinedata[0]), int(filelinedata[1]), filelinedata[2], int(filelinedata[3]),
                                     int(filelinedata[4]), int(filelinedata[5]), filelinedata[6], filelinedata[7],
                                     float(filelinedata[8]), int(filelinedata[9]), int(filelinedata[10]), float(filelinedata[11]), conrad_data)
                exomedepth_calls.append(edc)
    except IOError:
        print(f"Could not read ExomeDepth calls file {exomedepthfileloc}")
    finally:
        return exomedepth_calls


def read_exomedepth_collection(exomedepthfolder):
    """Reads a set of ExomeDepth call files

    Parameters
    ----------
    exomedepthfolder : str
        Path to folder containing ExomeDepth call files.

    Returns
    -------
    exomedepth_data : dict
        ExomeDepth CNV calls per sample.
    """
    ed_files = os.listdir(exomedepthfolder)
    ed_files = [edfile for edfile in ed_files if edfile.endswith(".csv")]

    exomedepth_data = {}
    for edfile in ed_files:
        sample_prefix = edfile.split("_")[0]
        if sample_prefix not in exomedepth_data:
            exomedepth_data[sample_prefix] = None
        sample_calls = read_exomedepth_calls(edfile)
        exomedepth_data[sample_prefix] = sample_calls
