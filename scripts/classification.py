import os
import argparse

# Import required classes
from classes.probe import Probe
from classes.exon import Exon
from classes.cnv import Cnv
from classes.arraycnv import ArrayCnv
from classes.conifercall import ConiferCall
from classes.exomedepthcall import ExomeDepthCall

# Import classification scripts
import classification.classification as clcl
#import classification.conifer as clfc
#import classification.exomedepth as clfed
#import classification.gatk as clfg

# Import parameters scripts
import parameters.parameters as parpar

# Import shared methods scripts
import shared_methods.shared_methods as smsm

# Import utils scripts
import utils.filereaders as ufr
import utils.filewriters as ufw


# Make parameter defining variables
TOOL_CHOICES = ["conifer", "exomedepth", "gatk"]
REQUIRED_PARAMS = {"conifer": ["arrayfile", "exonsfile", "infile", "probesfile", "samples"],
                   "exomedepth": ["arrayfile", "exonsfile", "infile", "probesfile", "samples"],
                   "gatk": ["arrayfile", "exonsfile", "indir", "probesfile", "samples"]}
OPTIONAL_PARAMS = {"conifer": ["numofexons", "numofprobes"],
                   "exomedepth": ["numofexons", "numofprobes"],
                   "gatk": ["numofexons", "numofprobes"]}
PARAM_TYPES = {"samples": "inputfile",
               "exonsfile": "inputfile",
               "probesfile": "inputfile",
               "arrayfile": "inputfile",
               "indir": "directory",
               "infile": "inputfile",
               "outprefix": "string",
               "numbofprobes": "integer",
               "numofexons": "integer",
               "percentoverlap": "integer",
               "cnvsize": "integer"}
TOOL_USAGE = {"conifer": "python classification.py -t conifer",
              "exomedepth": "python classification.py -t exomedepth",
              "gatk": "python classifiation.py -t gatk"}

# Make call translation tables for  the three CNV 
GATK_CALL_TRANSLATIONS = {'+': "CN Gain",
                          '-': "CN Loss"}
CONIFER_CALL_TRANSLATIONS = {"dup": "CN Gain",
                             "del": "CN Loss"}
EXOMEDEPTH_CALL_TRANSLATIONS = {"duplication": "CN Gain",
                                "deletion": "CN Loss"}


def main():
    """Performs the main work.

    Parameters
    ----------
    cmd_argvalues
        Command line parameter values
    """
    cmd_argvalues = parpar.get_classification_parameters(TOOL_CHOICES)
    incorrect_parameters = parpar.parameters_are_ok(cmd_argvalues, REQUIRED_PARAMS, PARAM_TYPES)

    if len(incorrect_parameters) == 0:
        print("...Reading the sample table...")
        sample_data = ufr.read_sample_table(cmd_argvalues["samples"])
        print("...Reading the probe file...")
        probe_data = ufr.read_probes_data(cmd_argvalues["probesfile"])
        print("...Reading the exon data...")
        exon_data = ufr.read_exon_data(cmd_argvalues["exonsfile"])
        print("...Reading the array data...")
        array_data = read_array_cnvs(cmd_argvalues["arrayfile"], exon_data)

        if len(sample_data) > 0 and len(probe_data) > 0 and len(exon_data) > 0 and len(array_data) > 0:
            # GATK4 CNV CLASSIFCATION
            if cmd_argvalues["tool"].upper() == "GATK":
                print("...Start combining GATK4 .called.seg files...")
                gatk4_cnv_data = gatk4_combine_seg_files(cmd_argvalues["indir"], sample_data, probe_data, exon_data, "")
                print(f"...Read GATK4 CNV data for {len(gatk4_cnv_data)} samples")
                print("...Evaluating GATK4 CNVs...")
                gatk4_cnv_data = gatk4_evaluate(sample_data, probe_data, exon_data, array_data, gatk4_cnv_data)
                #print(gatk4_cnv_data)
                #tmp_print_gatkcnvs(gatk4_cnv_data)
                print("...Classifying leftover array CNVs...")
                array_classify_leftovers(array_data)
                print("...Writing all GATK4 CNV classifications to file...")
                # write_cnv_results(gatk4_cnv_data, cmd_argvalues["output"]+"/gatk4_results.txt")
                write_cnv_results_2(gatk4_cnv_data, cmd_argvalues["output"], cmd_argvalues["filterneutrals"], cmd_argvalues["cnvsize"])
                # print("Writing Array CNVs with no GATK4 CNVs to file...")
                # leftoverpath = get_leftover_outpath(cmd_argvalues["output"])
                # write_array_leftovers(array_data, f"{leftoverpath}/array_leftovers.txt")

            # CONIFER CNV CLASSIFICATION
            if cmd_argvalues["tool"] == "conifer":
                print("...Start reading the Conifer data...")
                # conifer_data = ufr.read_conifer_calls(cmd_argvalues["indir"])
                conifer_data = read_conifer_data(cmd_argvalues["infile"], sample_data, exon_data, probe_data)
                print(f"...Read Conifer CNV data for {len(conifer_data)} samples...")

                print("...Evaluating Conifer CNVs...")
                # conifer_evaluate(conifer_data, array_data, cmd_argvalues["numofexons"], cmd_argvalues["numofprobes"], cmd_argvalues["percentoverlap"], CONIFER_CALL_TRANSLATIONS)
                evaluate_cnv_calls(conifer_data, array_data, cmd_argvalues["numofexons"], cmd_argvalues["numofprobes"], cmd_argvalues["percentoverlap"], CONIFER_CALL_TRANSLATIONS)

                print("...Classifying leftover array CNVs...")
                array_classify_leftovers(array_data)

                print("...Writing all Conifer classifications to file...")
                file_written = write_cnv_classifications(conifer_data, cmd_argvalues["output"], cmd_argvalues["filterneutrals"], cmd_argvalues["cnvsize"], CONIFER_CALL_TRANSLATIONS, "Conifer")
                print(f"...Output file written?: {file_written}...")

            # EXOMEDEPTH CNV CLASSIFICATION
            if cmd_argvalues["tool"] == "exomedepth":
                print("...Start reading the ExomeDepth data...")
                exomedepth_data = read_combined_exomedepth_data(cmd_argvalues["infile"], sample_data, exon_data, probe_data)
                print(f"...Read ExomeDepth CNV data for {len(exomedepth_data)} samples...")

                print("...Evaluating ExomeDepth CNVs...")
                evaluate_cnv_calls(exomedepth_data, array_data, cmd_argvalues["numofexons"], cmd_argvalues["numofprobes"], cmd_argvalues["percentoverlap"], EXOMEDEPTH_CALL_TRANSLATIONS)

                print("...Classifying leftover array CNVs...")
                array_classify_leftovers(array_data)

                print("...Writing all ExomeDepth classifications to file...")
                file_written = write_cnv_classifications(exomedepth_data, cmd_argvalues["output"], cmd_argvalues["filterneutrals"], cmd_argvalues["cnvsize"], EXOMEDEPTH_CALL_TRANSLATIONS, "ExomeDepth")
                print(f"...Outout file written?: {file_written}...")
            print("DONE!")


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


# ===================================================================================================================================================
def read_conifer_data(coniferfileloc, sampledata, exondata, probedata):
    """Read and reutrn conifer CNV call data.

    Parameters
    ----------
    coniferfileloc : str
        Path to file with Conifer CNV calls
    sampledata : dict
        Sample translation table
    exondata : dict
        Exon data per chromosome
    probedata : dict
        Probe data per chromosome

    Returns
    -------
    conifer_data : dict
        Conifer CNV calls with overlapping probes and exons per sample
    """
    conifer_data = {}
    try:
        with open(coniferfileloc, 'r') as coniferfile:
            next(coniferfile)
            for fileline in coniferfile:
                filelinedata = fileline.strip().split("\t")
                samplename = filelinedata[0].split(".")[0]

                # Get the pseudo sample name
                sample_pseudo = ""
                if samplename in sampledata:
                    sample_pseudo = sampledata[samplename]

                # Add the pseudo sample to the conifer_data
                if sample_pseudo != "" and sample_pseudo not in conifer_data:
                    conifer_data[sample_pseudo] = []

                # Create the ConiferCall, add overlapping exons and probes and add it to the data dict
                cfcall = ConiferCall(sample_pseudo, samplename, filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), filelinedata[4])
                cfcall.exons = clcl.add_exons(cfcall, exondata)
                cfcall.probes = clcl.add_probes(cfcall, probedata)
                conifer_data[sample_pseudo].append(cfcall)
    except IOError:
        print(f"Could not read the Conifer data at {coniferfileloc}")
    finally:
        return conifer_data


def evaluate_cnv_calls(tool_cnvs, array_cnvs, minexons, minprobes, minoverlap, calltranslationtable):
    for samplepseudo in tool_cnvs:
        for toolcnv in tool_cnvs[samplepseudo]:
            if toolcnv.cnv_sample_pseudo in array_cnvs:
                arraysamplecnvs = array_cnvs[toolcnv.cnv_sample_pseudo]
                classify_cnv(toolcnv, arraysamplecnvs, minexons, minprobes, minoverlap, calltranslationtable)


def conifer_evaluate(conifer_cnvs, array_cnvs, minexons, minprobes, minoverlap, calltranslationtable):
    for samplepseudo in conifer_cnvs:
        for cfcnv in conifer_cnvs[samplepseudo]:
            if cfcnv.cnv_sample_pseudo in array_cnvs:
                arraysamplecnvs = array_cnvs[cfcnv.cnv_sample_pseudo]
                classify_cnv(cfcnv, arraysamplecnvs, minexons, minprobes, minoverlap, calltranslationtable)
    return conifer_cnvs


def classify_cnv(cnvcall, arraycnvs, minexons, minprobes, minoverlap, translationtable):
    """Classify a GATK4, Conifer or ExomeDepth CNV call using an array CNV call.

    Parameters
    ----------
    cnvcall : Cnv or ConiferCall or ExomeDepthCall
        GATK4, Conifer or ExomeDepth CNV call to classify
    arraycnvs
        Called array CNVs for the same sample
    minexons : int
        Minimum number of exons required for a CNV call 
    minprobes : int
        Minimum number of probes required for a CNV call to be considered Array Informative
    minoverlap : int
        Minimum percentage overlap required for a CNV call to be considered a True Positive
    """
    match_with_array = False
    for arraycnv in arraycnvs:
        if cnvcall.cnv_chrom == arraycnv.cnv_chrom and cnvcall.cnv_overlap(arraycnv):
            cnv_overlap = cnvcall.get_percent_overlap(arraycnv.cnv_start, arraycnv.cnv_end)
            if cnv_overlap >= minoverlap:
                match_with_array = True
                cnvcall.classification = ["TRUE POSITIVE", f"{len(cnvcall.exons)}", f"{len(cnvcall.probes)}"]
                cnvcall.call_result = determine_call_result(cnvcall, arraycnv, translationtable)
                cnvcall.array_cnv = arraycnv
                cnvcall_hangovers = determine_hangover(cnvcall.cnv_start, cnvcall.cnv_end, arraycnv.cnv_start, arraycnv.cnv_end)
                cnvcall.left_hangover = cnvcall_hangovers[0]
                cnvcall.right_hangover = cnvcall_hangovers[1]
                arraycnv.wes_cnvs.append(cnvcall)

        # Determine the classification of the CNV call if there is no or too little overlap.
        if not match_with_array:
            if len(cnvcall.exons) >= minexons and len(cnvcall.probes) < minprobes:
                cnvcall.classification = ["ARRAY NON-INFORMATIVE", f"{len(cnvcall.exons)}", f"{len(cnvcall.probes)}"]
                cnvcall.call_result = "No array"
            if len(cnvcall.exons) < minexons and len(cnvcall.probes) >= minprobes:
                cnvcall.classification = ["WES NON-INFORMATIVE", f"{len(cnvcall.exons)}", f"{len(cnvcall.probes)}"]
                cnvcall.call_result = "No array"
            if len(cnvcall.exons) < minexons and len(cnvcall.probes) < minprobes:
                cnvcall.classification = ["ARRAY & WES NON-INFORMATIVE", f"{len(cnvcall.exons)}", f"{len(cnvcall.probes)}"]
                cnvcall.call_result = "No array"
            if len(cnvcall.exons) >= minexons and len(cnvcall.probes) >= minprobes:
                cnvcall.classification = ["FALSE POSITIVE", f"{len(cnvcall.exons)}", f"{len(cnvcall.probes)}"]
                cnvcall.call_result = "No array"


def determine_call_result(cnvcall, arraycnv, translationtable):
    """Determine and return the CNV call result.

    Parameters
    ----------
    gatkcnv : Cnv
        CNV called by GATK4
    arraycnv: ArrayCnv
        CNV called by the array
    translationtabel : dict
        Dict containing WES to array call translations

    Returns
    -------
    cnvcallresult : str
        CNV call result (Concordant, Discordant, Conflicting)
    """
    cnvcallresult = ""
    if cnvcall.cnv_call in translationtable:
        if translationtable[cnvcall.cnv_call] == arraycnv.cnv_call:
            cnvcallresult = "Concordant"
        else:
            cnvcallresult = "Conflicting"
    else:
        if cnvcall.cnv_call == '0':
            cnvcallresult = "Discordant"
    return cnvcallresult


# ===================================================================================================================================================
def read_exomedepth_sample_table(sample_ed_fileloc):
    """Read and return Sample to ExomeDepth table. Sample names should be the full sample names.

    Parameters
    ----------
    sample_ed_fileloc : str
        Path to Sample to ExomeDepth calls files table
    """
    sample_ed_table = {}
    try:
        with open(sample_ed_fileloc, 'r') as edsamplefile:
            next(edsamplefile)
            for fileline in edsamplefile:
                filelinedata = fileline.strip().split("\t")
                sample_ed_table[filelinedata[0]] = filelinedata[1]
    except IOError:
        print("Could not read sample to exomedepth file table")
    finally:
        return sample_ed_table


def combine_exomedepth_files(edfilesdir, sample_to_edfile, sampledata, probedata, exondata):
    indirfiles = os.listdir(edfilesdir)
    edfiles = [f"{edfilesdir}/{ifile}" for ifile in indirfiles if ifile.endswith(".csv")]

    combined_ed_file = {}
    for sample_name in sample_to_edfile:
        if sampledata[sample_name] not in combined_ed_file:
            combined_ed_file[sampledata[sample_name]] = []
        if sample_to_edfile[sample_name] in edfiles:
            combined_ed_file = read_exomedepth_data(sample_to_edfile[sample_name], combined_ed_file, sample_name, sampledata[sample_name], probedata, exondata)
    return combined_seg_file


def read_exomedepth_data(edfileloc, exomedepthdata, samplename, samplepseudo, probedata, exondata):
    try:
        with open(edfileloc, 'r') as edfile:
            next(edfile)
            for edline in edfile:
                edlinedata = edline.replace("\"", "").strip.split(",")

                edcall = ExomeDepthCall(samplename, samplepseudo, int(edlinedata[0]), int(edlinedata[1]), edlinedata[2], int(edlinedata[3]), int(edlinedata[4]), int(edlinedata[5]), edlinedata[6], edlinedata[7], float(edlinedata[8]), int(edlinedata[9]), int(edlinedata[10]), float(edlinedata[11]), edlinedata[12:])
                if samplename not in exomedepthdata:
                    exomedepthdata[samplename] = []
                exomedepthdata[samplename].append(edcall)
    except IOError:
        print(f"Could not read ExomeDepth calls file {edfileloc}")
    finally:
        return exomedepthdata


def read_combined_exomedepth_data(exdresultsfile, sampledata, exondata, probedata):
    """Read the combined ExomeDepth results file.

    Parameters
    ----------
    exd_resultsfile : str
        Path to combined ExomeDepth file
    sampledata : dict
        Sample translation table
    exondata : dict
        Exon data per chromosome
    probedata : dict
        Probe data per chromosome
    """
    exomedepthdata = {}
    try:
        with open(exdresultsfile, 'r') as exdfile:
            next(exdfile)
            for exdline in exdfile:
                exdlinedata = exdline.strip().split("\t")

                # Get the sample pseudo name
                sample_pseudo = ""
                if exdlinedata[0] in sampledata:
                    sample_pseudo = sampledata[exdlinedata[0]]

                # Add an entry in the exomedepth data
                if sample_pseudo != "" and sample_pseudo not in exomedepthdata:
                    exomedepthdata[sample_pseudo] = []

                # Create the ExomeDepthCall
                exdcall = ExomeDepthCall(exdlinedata[0], sample_pseudo, int(exdlinedata[1]), int(exdlinedata[2]), exdlinedata[3], int(exdlinedata[4]), int(exdlinedata[5]), int(exdlinedata[6]), f"chr{exdlinedata[7]}", exdlinedata[8], float(exdlinedata[9]), int(exdlinedata[10]), int(exdlinedata[11]), float(exdlinedata[12]), exdlinedata[13])
                exdcall.exons = clcl.add_exons(exdcall, exondata)
                exdcall.probes = clcl.add_probes(exdcall, probedata)
                exomedepthdata[sample_pseudo].append(exdcall)
    except IOError:
        print(f"Could not read combined ExomeDepth results file {exdresultsfile}")
    finally:
        return exomedepthdata


def write_cnv_classifications(toolcnvdata, outfileloc, filterneutrals, mincnvsize, calltranslationtable, toolname):
    """Write the classified CNV calls of a tool to a specified output file.

    Parameters
    ----------
    toolcnvdata : dict
        Classified CNVs for a specified tool
    outfileloc : str
        Path to write output file to
    filterneutrals : bool
        Whether to filter out neutral CNV calls (for GATK4 results)
    mincnvsize : int
        Minimal required CNV size
    calltranslationtable : dict
        Table to translate tool CNV calls to array CNV calls
    toolname : str
        Name of the tool used to call the classififed CNVs

    Returns
    -------
    wrote_file : bool
        True if output file has been written, False if not
    """
    wrote_file = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(f"Sample\t{toolname}_CNV\t{toolname}_Call\t{toolname}_Size\tArray_CNV\tArray_Call\tArray_Size\tHangover_L\tHangover_R\tCall_Result\tClassification\t#_Exons\t#_Probes\t{toolname}_genes\tArray_genes\t{toolname}_UGenes\tArray_UGenes\n")
            for samplename in toolcnvdata:
                for toolcnv in toolcnvdata[samplename]:
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
                    tool_genes = toolcnv.get_gene_names()
                    tool_gene_names = ":".join(tool_genes)
                    tool_ugene_names = tool_gene_names

                    # Replace the default values if there is an array CNV call for the tool CNV.
                    if toolcnv.array_cnv:
                        array_region = toolcnv.array_cnv.get_region()
                        array_call = toolcnv.array_cnv.cnv_call
                        array_size = toolcnv.array_cnv.get_length()
                        array_genes = toolcnv.array_cnv.get_gene_names()
                        array_gene_names = ":".join(array_genes)
                        left_hangover = toolcnv.left_hangover
                        right_hangover = toolcnv.right_hangover
                    
                    # Determine unique gatk and array gene names
                    if array_genes != "NA":
                        tool_ugene_names = ":".join(determine_unique_genes(tool_genes, array_genes))
                        array_ugene_names = ":".join(determine_unique_genes(array_genes, tool_genes))

                    # Check whether to filter out neutral GATK4 CNVs
                    if filterneutrals and toolcnv.cnv_call == '0':
                        write_line = False

                    # Check whether to filter out GATK4 CNVs smaller than a certain size
                    if mincnvsize is not None:
                        if toolcnv.get_length() < mincnvsize:
                            write_line = False

                    if write_line:
                        outfile.write(f"{samplename}\t{toolcnv.get_region()}\t{calltranslationtable[toolcnv.cnv_call]}\t{toolcnv.get_length()}\t"
                                      f"{array_region}\t{array_call}\t{array_size}\t{left_hangover}\t{right_hangover}\t"
                                      f"{toolcnv.call_result}\t{toolcnv.classification[0]}\t{toolcnv.classification[1]}\t"
                                      f"{toolcnv.classification[2]}\t{tool_gene_names}\t{array_gene_names}\t{tool_ugene_names}\t{array_ugene_names}\n")
        wrote_file = True
    except IOError:
        print(f"Could not write {toolname} CNV classification results to {outfileloc}")
    finally:
        return wrote_file


if __name__ == "__main__":
    main()
