from classes.conradcnv import ConradCnv
from classes.exon import Exon
from classes.gatkcall import GatkCall


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


def read_conrad_data(infileloc):
    """Read and return the Conrad CNVs. Conrad CNVs are saved per chromosome.

    Parameters
    ----------
    infileloc : str
        Path to file with Conrad CNV data
    """
    conrad_cnvs = {}
    list_index = 0
    try:
        with open(infileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                conradchrom = f"chr{filelinedata[0]}"

                if conradchrom not in conrad_cnvs:
                    conrad_cnvs[conradchrom] = []
                    list_index = 0
                conrad_cnv = ConradCnv(conradchrom, int(filelinedata[1]), int(filelinedata[2]), int(filelinedata[3]), filelinedata[4], filelinedata[5], list_index)

                if list_index > 0:
                    prev_index = list_index - 1
                    conrad_cnv.set_prev(conrad_cnvs[conradchrom][prev_index])
                    conrad_cnvs[conrad_chrom][prev_index].set_next(conrad_cnv)
                conrad_cnvs[conradchrom].append(conrad_cnv)
    except IOError:
        print(f"Could not read Conrad CNV file {infileloc}")
    finally:
        return conrad_cnvs


def add_qxte_exon_data(exonfileloc, conrad_cnv_data):
    """Add overlapping QXTE exons to Conrad CNVs.

    Parameters
    ----------
    exonfileloc : str
        Path to QXTE exon BED file
    conrad_cnv_data : dict
        Conrad CNV data
    """
    try:
        with open(exonfileloc, 'r') as exonfile:
            for exonline in exonfile:
                exonlinedata = exonline.strip().split("\t")

                # Form the exon chromosome name
                exonchromname = exonlinedata[0]
                if not exonchromname.startswith("chr"):
                    exonchromname = f"chr{exonchromname}"

                qxte_exon = Exon(exonchromname, int(exonlinedata[1]), int(exonlinedata[2]), exonlinedata[3])
                if exonchromname in conrad_cnv_data:
                    add_exons_to_conradcnv(conrad_cnv_data[exonchromname], qxte_exon)
    except IOError:
        print(f"Could not read {exonfileloc}")


def add_conrad_exon_data(exonfileloc, conrad_cnv_data):
    """Add Conrad exons to already read Conrad CNVs.

    Parameters
    ----------
    exonfileloc : str
        Path to file with Conrad exons
    conrad_cnv_data : dict
        Saved Conrad CNVs per chromosome
    """
    try:
        with open(exonfileloc, 'r') as cexonfile:
            next(cexonfile)
            for exonline in cexonfile:
                exonlinedata = exonline.strip().split("\t")

                conradchrom = f"chr{exonlinedata[0]}"
                cexon_start = int(exonlinedata[1])
                cexon_end = int(exonlinedata[2])
                cexon_names = exonlinedata[3].split(",")

                exon_names = []
                for exonname in cexon_names:
                    exon_names.append(exonname.split("_")[0])
                conrad_exon = ConradExon(conradchrom, cexon_start, cexon_end, exon_names)

                if conradchrom in conrad_cnv_data:
                    add_exons_to_conradcnv(conrad_cnv_data[conradchrom], conrad_exon)
    except IOError:
        print(f"Could not read Conrad exon file {exonfileloc}")
    finally:
        return conrad_cnv_data


def add_exons_to_conradcnv(conradchromcnvs, conradexon):
    """Add a specified Conrad exon to an overlapping Conrad CNV.

    Parameters
    ----------
    conradchromcnvs : list of ConradCnv
        Conrad CNVs located on a single chromosome
    conradexon : ConradExon
        Conrad exon to add to an overlapping Conrad CNV
    """
    for conradcnv in conradchromcnvs:
        if conrad_cnv_exon_overlap(conradcnv.conrad_start, conradcnv.conrad_end, conradexon.exon_start, conradexon.exon_end):
            conradcnv.add_conrad_exon(conradexon)
            break


def conrad_cnv_exon_overlap(c_cnvstart, c_cnvend, c_exonstart, c_exonend):
    """Check and return whether a specified Conrad exon overlaps with a specified Conrad CNV.

    Parameters
    ----------
    c_cnvstart : int
        Leftmost genomic position of the Conrad CNV
    c_cnvend : int
        Rightmost genomic position of the Conrad CNV
    c_exonstart : int
        Leftmost genomic position of the Conrad exon
    c_exonend : int
        Rightmost genomic position of the Conrad exon
    
    Returns
    -------
    bool
        True if Conrad exon overlaps with Conrad CNV, False if not
    """
    return c_cnvstart <= c_exonend and c_exonstart <= c_cnvend


def determine_gatk_conrad_overlaps(gatkcnvs, conradcnvdata):
    gatkcnvs_with_conrads = {}
    for samplename in gatkcnvs:
        for gatkcnv in gatkcnvs[samplename]:
            get_conrad_in_gatk(gatkcnv, conradcnvdata)

            if gatkcnv.num_of_conrad_cnvs() > 0:
                if gatkcnv.samplename not in gatkcnvs_with_conrads:
                    gatkcnvs_with_conrads[gatkcnv.samplename] = []
                gatkcnvs_with_conrads[gatkcnv.samplename].append(gatkcnv)
    return gatkcnvs_with_conrads


def get_conrad_in_gatk(gatkcnv, conradcnvdata):
    """Add overlapping Conrad CNVs to a specified GATK4 CNV. 

    Parameters
    ----------
    gatkcnv : GatkCall
        GATK4 CNV to add overlapping Conrad CNVs to
    conradcnvdata : dict
        Conrad CNV data
    """
    if gatkcnv.chrom in conradcnvdata:
        for conradcnv in conradcnvdata[gatkcnv.chrom]:
            if conrad_gatk_overlap(gatkcnv, conradcnv):
                gatkcnv.add_conrad_cnv(conradcnv)


def conrad_gatk_overlap(gatkcnv, conradcnv):
    """Determine and return whether a specified Conrad CNV overlaps with a specified GATK4 CNV

    Parameters
    ----------
    gatkcnv : GatkCall
        GATK4 CNV to check overlap against
    conradcnv : ConradCnv
        Conrad CNV to check overlap for
    """
    return gatkcnv.startpos <= conradcnv.conrad_end and conradcnv.conrad_start <= gatkcnv.endpos


def determine_gatk_filtered_by_conrad(gatkconraddata):
    """Determine whether the genenames in 

    Parameters
    ----------
    gatkconraddata : dict
        GATK4 CNVs overlapping with one or more Conrad CNVs
    """
    conrad_filtered_gatk = {}
    for samplename in gatkconraddata:
        for gatkcall in gatkconraddata[samplename]:
            gatkgenenames = gatkcall.gene_names
            conradgenenames = gatkcall.get_conrad_cnv_genes()

            shared_genenames = list(set(gatkgenenames) & set(conradgenenames))
            ugatk_genenames = list(set(gatkgenenames) - set(conradgenenames))
            uconrad_genenames = list(set(conradgenenames) - set(gatkgenenames))

            if len(ugatk_genenames) == 0:
                if samplename not in conrad_filtered_gatk:
                    conrad_filtered_gatk[samplename] = []
                conrad_filtered_gatk[samplename].append(gatkcall.get_region_str())
    return conrad_filtered_gatk


def get_gatk_region(gatkregion):
    """Split and return a GATK4 CNV region into three separate fields.

    Parameters
    ----------
    gatkregion : str
        GATK4 CNV region as chr:start-end

    Returns
    -------
    list of str and int
        Genome region split into chromosome, start and stop
    """
    gatk_chrom = gatkregion.split(":")[0]
    gatk_start = int(gatkregion.split(":")[1].split("-")[0])
    gatk_end = int(gatkregion.split(":")[1].split("-")[1])
    return [gatk_chrom, gatk_start, gatk_end]


def write_conrad_filtered_gatkcalls(gatkcalldata, filteredgatkcnvs, outfileloc):
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tGATK4_CNV\tGATK4_Call\tGATK4_Size\tArray_CNV\tArray_Call\tArray_Size\tHangover_L\tHangover_R\tCall_Result\tClassification\t#_Exons\t#_Probes\tGATK4_genes\tArray_genes\tGATK4_UGenes\tArray_UGenes\n")
            for samplename in gatkcalldata:
                for gatkcnv in gatkcalldata[samplename]:
                    if samplename in filteredgatkcnvs:
                        if gatkcnv.get_region_str() not in filteredgatkcnvs[samplename]:
                            outfile.write(f"{gatkcnv.call_line}\n")
                    else:
                        outfile.write(f"{gatkcnv.call_line}\n")
        file_written = True
    except IOError:
        print(f"Could not write GATK CNV calls to {outfileloc}")
    finally:
        return file_written


def write_gatk_conrad_overlaps(gatkconraddata, outfileloc):
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tChrom\tG_Start\tC_Start\tG_End\tC_End\n")
            for samplename in gatkconraddata:
                for gatkcnv in gatkconraddata[samplename]:
                    for ccnv in gatkcnv.conrad_cnvs:
                        outfile.write(f"{samplename}\t{gatkcnv.chrom}\t{gatkcnv.startpos}\t{ccnv.conrad_start}\t{gatkcnv.endpos}\t{ccnv.conrad_end}\n")
    except IOError:
        print(f"Could not write to {outfileloc}")

