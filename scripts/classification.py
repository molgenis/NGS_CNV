import os
import argparse

from probe import Probe
from exon import Exon
from cnv import Cnv
from arraycnv import ArrayCnv


def main(cmd_argvalues):
    """Performs the main work.

    Parameters
    ----------
    cmd_argvalues
        Command line parameter values
    """
    print("Reading the sample table...")
    sample_data = read_sample_table(cmd_argvalues["samples"])
    print("Reading the probe file...")
    probe_data = read_probes_data(cmd_argvalues["probes"])
    print("Reading the exon data...")
    exon_data = read_exon_data(cmd_argvalues["exons"])
    print("Reading the array data...")
    array_data = read_array_cnvs(cmd_argvalues["array"], exon_data)

    if len(sample_data) > 0 and len(probe_data) > 0 and len(exon_data) > 0 and len(array_data) > 0:
        if cmd_argvalues["tool"].upper() == "GATK":
            print("Start combining GATK4 .called.seg files...")
            gatk4_cnv_data = gatk4_combine_seg_files(cmd_argvalues["indir"], sample_data, probe_data, exon_data, "")
            print(f"Read GATK4 CNV data for {len(gatk4_cnv_data)} samples")
            print("Evaluating GATK4 CNVs...")
            gatk4_cnv_data = gatk4_evaluate(sample_data, probe_data, exon_data, array_data, gatk4_cnv_data)
            #print(gatk4_cnv_data)
            #tmp_print_gatkcnvs(gatk4_cnv_data)
            print("Classifying leftover array CNVs...")
            array_classify_leftovers(array_data)
            print("Writing all GATK4 CNV classifications to file...")
            # write_cnv_results(gatk4_cnv_data, cmd_argvalues["output"]+"/gatk4_results.txt")
            write_cnv_results_2(gatk4_cnv_data, cmd_argvalues["output"], cmd_argvalues["filterneutrals"], cmd_argvalues["cnvsize"])
            # print("Writing Array CNVs with no GATK4 CNVs to file...")
            # leftoverpath = get_leftover_outpath(cmd_argvalues["output"])
            # write_array_leftovers(array_data, f"{leftoverpath}/array_leftovers.txt")
            print("DONE!")


def get_args():
    """Define, receive and return command line parameters values.

    Returns
    -------
    dict
        Command line argument values
    """
    cmd_args = argparse.ArgumentParser()
    cmd_args.add_argument("-t", "--tool", type=str, dest="tool", help="Name of the tool")
    cmd_args.add_argument("-s", "--samples", type=str, dest="samples", help="Path to samples file")
    cmd_args.add_argument("-e", "--exons", type=str, dest="exons", help="Path to exons file")
    cmd_args.add_argument("-p", "--probes", type=str, dest="probes", help="Path to probes file")
    cmd_args.add_argument("-i", "--indir", type=str, dest="indir", help="Path to input folder or file")
    cmd_args.add_argument("-a", "--array", type=str, dest="array", help="Path to array calls file")
    cmd_args.add_argument("-o", "--output", type=str, dest="output", help="Path to write output file to")
    cmd_args.add_argument("-op", "--outprefix", type=str, dest="outprefix", help="Prefix to use for output files")
    cmd_args.add_argument("-np", "--number-of-probes", type=int, dest="numofprobes", default=10, help="Required minimum number of probes to be Array Informative")
    cmd_args.add_argument("-ne", "--number-of-exons", type=int, dest="numofexons", default=3, help="Required minimum number of exons to be WES Informative")
    cmd_args.add_argument("-po", "--percent-overlap", type=int, dest="percentoverlap", default=50, help="Percentage overlap required between Array en WES CNVs")
    cmd_args.add_argument("-fn", "--filter-neutrals", dest="filterneutrals", action="store_true", help="Filter out neutrals?")
    cmd_args.add_argument("-cs", "--cnv-size", type=int, dest="cnvsize", help="Minimum GATK4 CNV size")
    return vars(cmd_args.parse_args())


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

                if exonlinedata[0] not in exon_data:
                    exon_data[exonlinedata[0]] = []
                exon_data[exonlinedata[0]].append(Exon(exonlinedata[0], int(exonlinedata[1]), int(exonlinedata[2]), exonlinedata[3]))
    except IOError:
        print(f"Could not open exon file {exonfileloc}")
    finally:
        return exon_data


def read_array_cnvs(arrayfileloc, exondata):
    """Read the array CNV data.

    Parameters
    ----------
    arrayfileloc : str
        Path to array CNV file
    exondata : dict
        Exon data

    Returns
    -------
    arraydata : dict
        Array CNV data
    """
    arraydata = {}
    try:
        with open(arrayfileloc, 'r') as arrayfile:
            next(arrayfile)
            for arrayline in arrayfile:
                arraylinedata = arrayline.strip().split("\t")
                cnv_region = array_get_cnv_region(arraylinedata[1])
                array_cnv = ArrayCnv(arraylinedata[0], cnv_region[0], cnv_region[1], cnv_region[2], arraylinedata[2], int(arraylinedata[3]), int(arraylinedata[4]), int(arraylinedata[5]), arraylinedata[6])
                array_cnv.exons = gatk4_add_exons(array_cnv, exondata)

                if arraylinedata[0] not in arraydata:
                    arraydata[arraylinedata[0]] = []
                arraydata[arraylinedata[0]].append(array_cnv)
    except IOError:
        print(f"Could not read array CNV file {arrayfileloc}")
    finally:
        return arraydata


def gatk4_combine_seg_files(segfilesdir, sampledata, probedata, exondata, outfilepath):
    """Combine seg files and output a new combined file.

    Parameters
    ----------
    segfilesdir : str
        Path to directory containing GATK4 CCRS .called.seg files
    outfilepath : str
        Path to write output file to

    Returns
    -------
    combined_seg_file : dict
        Combined GATK4 segment data
    """
    skiplines = ["@HD", "@SG", "CONTIG"]
    indirfiles = os.listdir(segfilesdir)
    calledsegfiles = [f"{segfilesdir}/{ifile}" for ifile in indirfiles if ifile.endswith(".called.seg")]

    combined_seg_file = {}
    for segfile in calledsegfiles:
        combined_seg_file = gatk4_read_calledseg_file(segfile, combined_seg_file, skiplines, sampledata, probedata, exondata)
    return combined_seg_file


def gatk4_read_calledseg_file(segfileloc, combinedsegdata, linestoskip, sampletable, probedata, exondata):
    """Read and add segment file to combined segment file data.

    Parameters
    ----------
    segfileloc : str
        Path to segment file to process
    combinedsegdata : dict
        Collected segment data from multiple segment files
    linetoskip : list of str
        Starts of lines to skip in segment files
    probedata : dict
        Probes per chromosome
    exondata : dict
        Exons per chromosome

    Returns
    -------
    combinedsegdata : dict
        Updated collected segment data from multiple segment files
    """
    try:
        with open(segfileloc, 'r') as segfile:
            segsample = ""
            sample_pseudo = ""
            for segline in segfile:
                seglinedata = segline.strip().split("\t")

                if seglinedata[0] not in linestoskip:
                    if seglinedata[0] == "@RG":
                        segsample = seglinedata[-1].split(":")[1]
                        sample_pseudo = get_pseudo_sample(segsample, sampletable)
                        if sample_pseudo not in combinedsegdata:
                            combinedsegdata[sample_pseudo] = []
                    else:
                        if sample_pseudo in combinedsegdata:
                            cnv_to_add = Cnv(segsample, sample_pseudo, seglinedata[0], int(seglinedata[1]), int(seglinedata[2]), float(seglinedata[3]), float(seglinedata[4]), seglinedata[5])
                            cnv_to_add.probes = gatk4_add_probes(cnv_to_add, probedata)
                            cnv_to_add.exons = gatk4_add_exons(cnv_to_add, exondata)
                            combinedsegdata[sample_pseudo].append(cnv_to_add)
    except IOError:
        print(f"Could not read GATK4 CCRS .called.seg file {segfileloc}")
    finally:
        return combinedsegdata


def get_pseudo_sample(samplename, sampletable):
    """Returns the pseudo sample name of the .

    Parameters
    ----------
    samplename : str
        Sample name to obtain pseudo name for
    sampletable : dict
        Sample table with sample name and pseudo names

    Returns
    -------
    str
        Pseudo sample name if available, empty name if not
    """
    if samplename in sampletable:
        return sampletable[samplename]
    return ""


def array_get_cnv_region(arrayregion):
    """Parse and return the region of an array CNV.

    Parameters
    ----------
    arrayregion : str
        

    Returns
    -------
    list of str and int
        Parsed CNV region data
    """
    regiondata = arrayregion.split(":")
    chrom = regiondata[0]
    startpos = int(regiondata[1].split("-")[0].replace(",", ""))
    endpos = int(regiondata[1].split("-")[1].replace(",", ""))
    return [chrom, startpos, endpos]


def gatk4_add_probes(cnv, probes):
    """Identify probes overlapping with the GATK CNV.

    Parameters
    ----------
    cnv : Cnv
        GATK4 CNV
    probes : dict
        Probe data
    """
    probes_to_add = []
    chromprobes = probes[cnv.cnv_chrom]
    for chromprobe in chromprobes:
        if cnv.segment_overlap(chromprobe.probe_start, chromprobe.probe_end):
            probes_to_add.append(chromprobe)
    return probes_to_add


def gatk4_add_exons(cnv, exons):
    exons_to_add = []
    chromexons = exons[cnv.cnv_chrom]
    for chromexon in chromexons:
        if cnv.segment_overlap(chromexon.exon_start, chromexon.exon_end):
            exons_to_add.append(chromexon)
    return exons_to_add


def gatk4_evaluate(sample_data, probe_data, exon_data, array_cnvs, gatk4_cnvs):
    """Evaluate GATK4 CNV calls with array calls per sample.

    Parameters
    ----------
    sample_data : dict
        Sample translation table
    probe_data : dict
        Probe locations per chromosome
    gatk4_cnvs : str
        Path to file with combined .called.seg files data
    """
    for samplename in gatk4_cnvs:
        print(f"Evaluating sample {samplename}")
        for gatk4cnv in gatk4_cnvs[samplename]:
            if gatk4cnv.cnv_sample_pseudo in array_cnvs:
                arraysamplecnvs = array_cnvs[gatk4cnv.cnv_sample_pseudo]
                gatk4_classify_cnv(gatk4cnv, arraysamplecnvs)
    return gatk4_cnvs


def gatk4_classify_cnv(gatkcnv, arraycnvs):
    """Classify and a CNV.

    Parameters
    ----------
    gatkcnv : Cnv
        GATK4 called CNV
    arraycnvs : list of ArrayCnv
        Sample Arrays CNV
    """
    match_with_array = False
    classification = ""
    callresult = ""
    for arraycnv in arraycnvs:
        if gatkcnv.cnv_chrom == arraycnv.cnv_chrom and gatkcnv.cnv_overlap(arraycnv):
            match_with_array = True
            gatkcnv.array_cnv = arraycnv
            gatk_hangovers = determine_hangover(gatkcnv.cnv_start, gatkcnv.cnv_end, arraycnv.cnv_start, arraycnv.cnv_end)
            gatkcnv.left_hangover = gatk_hangovers[0]
            gatkcnv.right_hangover = gatk_hangovers[1]
            arraycnv.wes_cnvs.append(gatkcnv)

            if gatkcnv.num_of_exons() >= 3 and arraycnv.number_of_probes >= 10:
                # classification = "POSITIVE"
                classification = ["POSITIVE", f"{gatkcnv.num_of_exons()}", f"{arraycnv.num_of_probes()}"]
                callresult = gatk4_determine_call_result(gatkcnv, arraycnv)
            elif gatkcnv.num_of_exons() >= 3 and arraycnv.number_of_probes < 10:
                # classification = f"POSITIVE ARRAY NON-INFORMATIVE ()"
                classification = ["POSITIVE ARRAY NON-INFORMATIVE", f"{gatkcnv.num_of_exons()}", f"{arraycnv.num_of_probes()}"]
                callresult = gatk4_determine_call_result(gatkcnv, arraycnv)
            elif gatkcnv.num_of_exons() < 3 and arraycnv.number_of_probes >= 10:
                # classification = f"POSITIVE WES NON-INFORMATIVE ()"
                classification = ["POSITIVE WES NON-INFORMATIVE", f"{gatkcnv.num_of_exons()}", f"{arraycnv.num_of_probes()}"]
                callresult = gatk4_determine_call_result(gatkcnv, arraycnv)
            elif gatkcnv.num_of_exons() < 3 and arraycnv.number_of_probes < 10:
                # classification = f"POSITIVE ARRAY & WES NON-INFORMATIVE ()"
                classification = ["POSITIVE ARRAY & WES NON-INFORMATIVE", f"{gatkcnv.num_of_exons()}", f"{arraycnv.num_of_probes()}"]
                callresult = gatk4_determine_call_result(gatkcnv, arraycnv)
    if not match_with_array:
        if gatkcnv.num_of_exons() >= 3 and gatkcnv.num_of_probes() < 10:
            # classification = f"ARRAY NON-INFORMATIVE ({gatkcnv.num_of_probes()} probes)"
            classification = ["ARRAY NON-INFORMATIVE", f"{gatkcnv.num_of_exons()}", f"{gatkcnv.num_of_probes()}"]
            callresult = "No array"
        elif gatkcnv.num_of_exons() < 3 and gatkcnv.num_of_probes() >= 10:
            # classification = f"WES NON-INFORMATIVE ({gatkcnv.num_of_exons()} exons)"
            classification = ["WES NON-INFORMATIVE", f"{gatkcnv.num_of_exons()}", f"{gatkcnv.num_of_probes()}"]
            callresult = "No array"
        elif gatkcnv.num_of_exons() < 3 and gatkcnv.num_of_probes() < 10:
            # classification = f"ARRAY & WES NON-INFORMATIVE ({} exons ; {} probes)"
            classification = ["ARRAY & WES NON-INFORMATIVE", f"{gatkcnv.num_of_exons()}", f"{gatkcnv.num_of_probes()}"]
            callresult = "No array"
        elif gatkcnv.num_of_exons() >= 3 and gatkcnv.num_of_probes() >= 10:
            # classification = f"FALSE POSITIVE ({gatkcnv.num_of_exons()} exons ; {gatkcnv.num_of_probes()} probes)"
            classification = ["FALSE POSITIVE", f"{gatkcnv.num_of_exons()}", f"{gatkcnv.num_of_probes()}"]
            callresult = "No array"
    gatkcnv.classification = classification
    gatkcnv.call_result = callresult


def array_classify_leftovers(arraydata):
    """Classify and return classifications for leftover array CNVs.

    Parameters
    ----------
    arraydata : dict
        Array CNVs

    Returns
    -------
    arraydata : dict
        Updated array data
    """
    for arraysample in arraydata:
        for arraycnv in arraydata[arraysample]:
            if arraycnv.num_of_wes_cnvs() < 0:
                arraycnv.classification = "FALSE NEGATIVE"
    return arraydata


def gatk4_determine_call_result(gatkcnv, arraycnv):
    """Determine and return the CNV call result.

    Parameters
    ----------
    gatkcnv : Cnv
        CNV called by GATK4
    arraycnv: ArrayCnv
        CNV called by the array

    Returns
    -------
    cnvcallresult : str
        CNV call result (Concordant, Discordant, Conflicting)
    """
    cnvcallresult = ""
    gatk_call_translate = {}
    gatk_call_translate['+'] = "CN Gain"
    gatk_call_translate['-'] = "CN Loss"
    
    if gatkcnv.cnv_call in gatk_call_translate:
        if gatk_call_translate[gatkcnv.cnv_call] == arraycnv.cnv_call:
            cnvcallresult = "Concordant"
        else:
            cnvcallresult = "Conflicting"
    else:
        if gatkcnv.cnv_call == '0':
            cnvcallresult = "Discordant"
    return cnvcallresult


def write_cnv_results(gatkcnvdata, outfileloc):
    """Write classified GATK4 CNVs to an output file.

    Parameters
    ----------
    gatkcnvdata : dict
        Classified GATK4 CNVs
    outfileloc : str
        Path to write output file to
    """
    gatk_call_translations = {}
    gatk_call_translations["+"] = "CN Gain"
    gatk_call_translations["-"] = "CN Loss"
    gatk_call_translations["0"] = "Neutral"

    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tGATK4_CNV\tGATK4_Call\tArray_CNV\tArray_Call\tHangover_L\tHangover_R\tCall_Result\tClassification\t#_Exons\t#_Probes\n")
            for samplename in gatkcnvdata:
                for gatkcnv in gatkcnvdata[samplename]:
                    array_region = "NA"
                    array_call = "NA"
                    left_hangover = "NA"
                    right_hangover = "NA"

                    if gatkcnv.array_cnv is not None:
                        array_region = gatkcnv.array_cnv.get_region()
                        array_call = gatkcnv.array_cnv.cnv_call
                        left_hangover = gatkcnv.left_hangover
                        right_hangover = gatkcnv.right_hangover

                    outfile.write(f"{samplename}\t{gatkcnv.get_region()}\t{gatk_call_translations[gatkcnv.cnv_call]}\t"
                                  f"{array_region}\t{array_call}\t{left_hangover}\t{right_hangover}\t{gatkcnv.call_result}\t"
                                  f"{gatkcnv.classification[0]}\t{gatkcnv.classification[1]}\t{gatkcnv.classification[2]}\n")
    except IOError:
        print(f"Could not write output file {outfileloc}")


def write_cnv_results_2(gatkcnvdata, outfileloc, filterneutrals, mincnvsize):
    gatk_call_translations = {}
    gatk_call_translations["+"] = "CN Gain"
    gatk_call_translations["-"] = "CN Loss"
    gatk_call_translations["0"] = "Neutral"

    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tGATK4_CNV\tGATK4_Call\tGATK4_Size\tArray_CNV\tArray_Call\tArray_Size\tHangover_L\tHangover_R\tCall_Result\tClassification\t#_Exons\t#_Probes\tGATK4_genes\tArray_genes\tGATK4_UGenes\tArray_UGenes\n")
            for samplename in gatkcnvdata:
                for gatkcnv in gatkcnvdata[samplename]:
                    write_line = True

                    # Make default values so writing to output file always works.
                    array_region = "NA"
                    array_call = "NA"
                    array_size = "NA"
                    array_genes = "NA"
                    array_gene_names = "NA"
                    array_ugene_names = "NA"
                    left_hangover = "NA"
                    right_hangover = "NA"
                    gatk_genes = gatkcnv.get_gene_names()
                    gatk_gene_names = ":".join(gatk_genes)
                    gatk_ugene_names = gatk_gene_names

                    # Replace the default values if there is an array CNV call for the GATK4 CNV.
                    if gatkcnv.array_cnv is not None:
                        array_region = gatkcnv.array_cnv.get_region()
                        array_call = gatkcnv.array_cnv.cnv_call
                        array_size = gatkcnv.array_cnv.get_length()
                        array_genes = gatkcnv.array_cnv.get_gene_names()
                        array_gene_names = ":".join(array_genes)
                        left_hangover = gatkcnv.left_hangover
                        right_hangover = gatkcnv.right_hangover

                    # Determine unique gatk and array gene names
                    if array_genes != "NA":
                        gatk_ugene_names = ":".join(determine_unique_genes(gatk_genes, array_genes))
                        array_ugene_names = ":".join(determine_unique_genes(array_genes, gatk_genes))

                    # Check whether to filter out neutral GATK4 CNVs
                    if filterneutrals and gatkcnv.cnv_call == '0':
                        write_line = False

                    # Check whether to filter out GATK4 CNVs smaller than a certain size
                    if mincnvsize is not None:
                        if gatkcnv.get_length() < mincnvsize:
                            write_line = False

                    # Check whether to write to the output file
                    if write_line:
                        outfile.write(f"{samplename}\t{gatkcnv.get_region()}\t{gatk_call_translations[gatkcnv.cnv_call]}\t{gatkcnv.get_length()}\t"
                                      f"{array_region}\t{array_call}\t{array_size}\t{left_hangover}\t{right_hangover}\t"
                                      f"{gatkcnv.call_result}\t{gatkcnv.classification[0]}\t{gatkcnv.classification[1]}\t"
                                      f"{gatkcnv.classification[2]}\t{gatk_gene_names}\t{array_gene_names}\t{gatk_ugene_names}\t{array_ugene_names}\n")
    except IOError:
        print(f"Could not write to output file: {outfileloc}")


def write_array_leftovers(arraycnvdata, outfileloc):
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tArray_CNV\tArray_Call\tClassification\n")
            
            for samplename in arraycnvdata:
                for arraycnv in arraycnvdata[samplename]:
                    if arraycnv.num_of_wes_cnvs() == 0 and arraycnv.classification == "FALSE NEGATIVE":
                        outfile.write(f"{samplename}\t{arraycnv.get_region()}\t{arraycnv.cnv_call}\t{arraycnv.classification}\n")
    except IOError:
        print(f"Could not write leftover array cnvs to {outfileloc}")


def determine_hangover(gatkcnvstart, gatkcnvend, arraycnvstart, arraycnvend):
    """Determine and return the left and right hangovers.

    Parameters
    ----------
    gatkcnvstart : int
        Leftmost genomic position of the GATK4 CNV
    gatkcnvend : int
        Rightmost genomic position of the GATK4 CNV
    arraycnvstart : int
        Leftmost genomic position of the array CNV
    arraycnvend : int
        Rightmost genomic position of the array CNV

    Returns
    -------
    hangovers : list of int
        Left and right hangover
    """
    hangovers = []
    hangovers.append(gatkcnvstart - arraycnvstart)
    hangovers.append(gatkcnvend - arraycnvend)
    return hangovers


def tmp_print_gatkcnvs(gatkcnvdata):
    for samplename in gatkcnvdata:
        print(f"CNVs for sample {samplename}")
        for gatkcnv in gatkcnvdata[samplename]:
            print(f"CNV {gatkcnv.get_region()} classified as {gatkcnv.classification[0]} ({gatkcnv.classification[1]} exons ; {gatkcnv.classification[2]} probes)")


def get_leftover_outpath(outfilepath):
    tmp_folderpath = outfilepath.split("/")
    if len(tmp_folderpath) > 1:
        return "/".join(tmp_folderpath[0:-1])
    return ""


def determine_unique_genes(genelist1, genelist2):
    if genelist1 is None:
        if genelist2 is None:
            return []
        else:
            return genelist2
    elif genelist2 is None:
        if genelist1 is None:
            return []
        else:
            return genelist1
    return [genename for genename in genelist1 if genename not in genelist2]


if __name__ == "__main__":
    cli_args = get_args()
    main(cli_args)
