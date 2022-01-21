import argparse
from ccrscall import CcrsCall
from read_combined_ccrs import read_combined_ccrs
from solverdcall import SolveRdCall
from solverdexomedepthcall import SolveRdExomeDepthCall
from erngene import ErnGene


CONIFER_TO_GATK4 = {"dup": '+',
                    "del": '-'}
EXOMEDEPTH_TO_GATK4 = {"duplication": '+',
                       "deletion": '-'}
CLINCNV_TO_GATK4 = {"DUP": '+',
                    "DEL": '-'}
VARGENIUS_TO_GATK = {"DUP": '+',
                     "DEL": '-'}
GATK4_CALLS_TRANSLATION = {'+': ["dup", "duplication", "DUP"],
                           '-': ["del", "deletion", "DEL"]}


def get_params():
    """Define, receive and return set CLI parameters."""
    TOOL_CHOICES = ["clincnv", "conifer", "exomedepth", "vargenius"]
    TYPE_CHOICES = ["overlap", "genes", "both"]
    compare_args = argparse.ArgumentParser()
    compare_args.add_argument("-g", "--gatk4-calls", type=str, required=True, dest="gatk4-calls", help="Path to combined CCRS file")
    compare_args.add_argument("-c", "--calls", required=True, dest="calls", help="Path to file with calls to compare to")
    compare_args.add_argument("-ct", "--callstool", required=True, choices=TOOL_CHOICES, dest="callstool", help="Tool used to make the other calls")
    compare_args.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", help="Path to write output to")
    compare_args.add_argument("-op", "--output-prefix", type=str, required=True, dest="output-prefix", help="Output prefix to use for")
    compare_args.add_argument("-p", "--percentage-overlap", type=float, dest="percentage-overlap", default=80.0, help="Minimum required percentage overlap")
    compare_args.add_argument("-e", "--ern-gene-list", type=str, dest="ern-gene-list", help="Path to ERN gene list")
    compare_args.add_argument("-t", "--type-of-comparison", type=str, choices=TYPE_CHOICES, dest="type-of-comparison", help="")
    return vars(compare_args.parse_args())


def read_ern_gene_list(ernfileloc):
    """Read the ERN gene list file."""
    ern_data = {}
    try:
        with open(ernfileloc, 'r') as ernfile:
            for fileline in ernfile:
                filelinedata = fileline.strip().split("\t")

                # Add the ERN gene to the data
                if filelinedata[1] not in ern_data:
                    ern_data[filelinedata[1]] = []
                ern_data[filelinedata[1]].append(ErnGene(filelinedata[0], filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), int(filelinedata[4]), int(filelinedata[5]), int(filelinedata[6]), int(filelinedata[8])))
    except IOError:
        print("Could not read ERN file")
    finally:
        return ern_data


def add_ern_genes_to_cnvcalls(ccrsdata, erndata):
    """Determine which ERN genes overlap with CCRS calls."""
    for samplename in ccrsdata:
        for chromname in ccrsdata[samplename]:
            for ccrscall in ccrsdata[samplename][chromname]:
                get_ern_genes_for_ccrscall(ccrscall, erndata[ccrscall.ccrs_chrom])


def get_ern_genes_for_ccrscall(ccrscall, ernchromdata):
    """Get the overlapping ERN genes for a single CCRS call."""
    for erngene in ernchromdata:
        if ccrs_erngene_overlap(ccrscall, erngene):
            ccrscall.ern_genes.append(erngene)


def ccrs_erngene_overlap(ccrscall, erngene):
    """Determine whether there is an overlap between the CCRS call and Conrad CNV."""
    if ccrscall.ccrs_chrom == erngene.ern_chrom:
        return ccrscall.ccrs_start <= erngene.ern_padded_stop and erngene.ern_padded_start <= ccrscall.ccrs_end
    return False


def read_exomedepth(exomedepthfileloc):
    """Read Solve-RD ExomeDepth file.

    Parameters
    ----------
    exomedepthfileloc : str
        Path to the Solve-RD ExomeDepth file

    Returns
    -------
    exomedepth_calls : dict
        Solve-RD 
    """
    exomedepth_calls = {}
    try:
        with open(exomedepthfileloc, 'r') as exdfile:
            next(exdfile)
            for fileline in exdfile:
                filelinedata = fileline.strip().split("\t")

                # Add the ExomeDepth call
                if filelinedata[0] not in exomedepth_calls:
                    exomedepth_calls[filelinedata[0]] = {}
                if filelinedata[9] not in exomedepth_calls[filelinedata[0]]:
                    exomedepth_calls[filelinedata[0]][filelinedata[9]] = []
                exomedepth_calls[filelinedata[0]][filelinedata[9]].append(SolveRdExomeDepthCall(filelinedata[0], filelinedata[9], int(filelinedata[10]), int(filelinedata[11]), int(filelinedata[27]), filelinedata[13], int(filelinedata[14])))
    except IOError:
        print("Could not read ExomeDepth calls file")
    finally:
        return exomedepth_calls


def read_conifer(coniferfileloc):
    """Read Solve-RD Conifer file.

    Parameters
    ----------
    coniferfileloc : str
        Path to the Solve-RD Conifer file

    Returns
    -------
    conifer_calls : dict
        Solve-RD Conifer calls per sample per chromosome
    """
    conifer_calls = {}
    try:
        with open(coniferfileloc, 'r') as conffile:
            next(conffile)
            for fileline in conffile:
                filelinedata = fileline.strip().split("\t")
                
                if filelinedata[1] not in conifer_calls:
                    conifer_calls[filelinedata[1]] = {}
                if filelinedata[8] not in conifer_calls[filelinedata[1]]:
                    conifer_calls[filelinedata[1]][filelinedata[8]] = []
                conifer_calls[filelinedata[1]][filelinedata[8]].append(SolveRdCall(filelinedata[1], filelinedata[8], int(filelinedata[9]), int(filelinedata[10]), int(filelinedata[11]), filelinedata[12]))
    except IOError:
        print("Could not read Conifer file")
    finally:
        return conifer_calls


def read_clincnv(clincnvfileloc):
    """Read the Solve-RD ClinCNV file.

    Parameters
    ----------
    clincnvfileloc : str
        Path to the Solve-RD ClinCNV file
    """
    clincnv_calls = {}
    try:
        with open(clincnvfileloc, 'r') as clincnvfile:
            next(clincnvfile)
            for fileline in clincnvfile:
                filelinedata = fileline.strip().split("\t")

                if filelinedata[1] not in clincnv_calls:
                    clincnv_calls[filelinedata[1]] = {}
                if filelinedata[10] not in clincnv_calls[filelinedata[1]]:
                    clincnv_calls[filelinedata[1]][filelinedata[10]] = []
                clincnv_calls[filelinedata[1]][filelinedata[10]].append(SolveRdCall(filelinedata[1], filelinedata[10], int(filelinedata[11]), int(filelinedata[12]), int(filelinedata[13]), filelinedata[24]))
    except IOError:
        print("Could not read ClinCNV file")
    finally:
        return clincnv_calls


def read_vargenius(vargeniusfileloc):
    """Read the Solve-RD VarGenius file.

    Parameters
    ----------
    vargeniusfileloc : str
        
    """
    vargenius_calls = {}
    try:
        with open(vargeniusfileloc, 'r') as vargeniusfile:
            next(vargeniusfile)
            for fileline in vargeniusfile:
                filelinedata = fileline.strip().split("\t")

                if filelinedata[1] not in vargenius_calls:
                    vargenius_calls[filelinedata[1]] = {}
                if filelinedata[14] not in vargenius_calls[filelinedata[1]]:
                    vargenius_calls[filelinedata[1]][filelinedata[14]] = []
                vargenius_calls[filelinedata[1]][filelinedata[14]].append(SolveRdCall(filelinedata[1], filelinedata[14], int(filelinedata[15]), int(filelinedata[16]), int(filelinedata[19]), filelinedata[17]))
    except IOError:
        print("Could not read VareGenius file")
    finally:
        return vargenius_calls


def compare_gatk_solverd_overlap(gatk4calls, exomedepthcalls, minpercoverlap, calltranslationtable):
    """Compare GATK4 and ExomeDepth CNV calls per sample."""
    shared_samples = set(exomedepthcalls.keys()) & set(gatk4calls.keys())
    for samplename in shared_samples:
        for chromname in exomedepthcalls[samplename]:
            if chromname in gatk4calls[samplename]:
                for exdcall in exomedepthcalls[samplename][chromname]:
                    determine_overlapping_ccrscalls(exdcall, gatk4calls[samplename][chromname], minpercoverlap, calltranslationtable)


def determine_overlapping_ccrscalls(othercall, chrom_ccrscalls, minpercoverlap, calltranslationtable):
    """Determine the CCRS calls overlapping with a single other Solve-RD call."""
    for ccrscall in chrom_ccrscalls:
        if determine_gatk_other_overlap(ccrscall.ccrs_start, ccrscall.ccrs_end, othercall.cnv_start, othercall.cnv_end) >= minpercoverlap:
            if calltranslationtable[othercall.cnv_call] == ccrscall.ccrs_call:
                othercall.overlapping_ccrs.append(ccrscall)
            else:
                othercall.overlapping_ccrs_2.append(ccrscall)


def determine_gatk_other_overlap(gatkstart, gatkend, otherstart, otherend):
    """Determine the percentage overlap between a group call and the group representative.

    Parameters
    ----------
    gatkstart : int
        Starting position of the GATK4 call
    gatkend : int
        Ending position of the GATK4 call
    otherstart : int
        Starting position of the other Solve-RD tool call
    otherend : int
        Ending position of the other Solve-RD tool call
    """
    gatk_length = gatkend - gatkstart
    overlapnum = max(0, min(otherend, gatkend) - max(otherstart, gatkstart))
    return round((overlapnum/gatk_length)*100, 2)


def main():
    compare_params = get_params()
    compare_params["outdir"] = compare_params["outdir"]+"/" if not compare_params["outdir"].endswith("/") else compare_params["outdir"]
    outpath = compare_params["outdir"] + compare_params["output-prefix"]

    # Read the CCRS calls
    print("[-READING CCRS ERN CALLS-]")
    ccrs_data = read_combined_ccrs(compare_params["gatk4-calls"])
    ccrs_header = ccrs_data[0]
    ccrs_calls = ccrs_data[1]

    # Perform the comparison with the ExomeDepth calls
    if compare_params["callstool"] == "exomedepth":
        print("[-READING SOLVE-RD EXOMEDEPTH CALLS-]")
        exd_calls = read_exomedepth(compare_params["calls"])

        print("[-COMPARING GATK4 AND SOLVE-RD EXOMEDEPTH CALLS-]")
        compare_gatk_solverd_overlap(ccrs_calls, exd_calls, compare_params["percentage-overlap"], EXOMEDEPTH_TO_GATK4)
        unique_gatk_samples = set(ccrs_calls.keys()) - set(exd_calls.keys())
        unique_exd_samples = set(exd_calls.keys()) - set(ccrs_calls.keys())

        # Write the overlaps to an output file
        print("[-WRITING EXOMEDEPTH CALLS WITH WHICH CCRS CALLS OVERLAP TO FILE-]")
        wrote_file = write_overlapping_calls(f"{outpath}_overlapping.txt", exd_calls, "ExomeDepth")
        print(f"...Wrote overlapping calls to file?: {wrote_file}...")

        # Write the unique sample names to output files
        print("[-WRITING UNIQUE GATK AND EXOMEDEPTH SAMPLE NAMES TO FILE-]")
        wrote_gatk4 = write_unique_samples(f"{outpath}_unique_gatk_samples.txt", unique_gatk_samples)
        print(f"...Wrote unique GATK4 sample names?: {wrote_gatk4}...")
        wrote_exd = write_unique_samples(f"{outpath}_unique_exomedepth_samples.txt", unique_exd_samples)
        print(f"...Wrote unique ExomeDepth samples names?: {wrote_exd}...")

        # Display the numbers
        count_samples_calls_without_ccrs(exd_calls)

        # Check if calls should also be compared based on ERN genes
        #if compare_params["type-of-comparison"] == "genes" or compare_params["type-of-comparison"] == "both":
        #    print("[-ADDING ERN GENES TO GATK4 CALLS-]")
        #    add_ern_genes_to_cnvcalls(ccrs_calls, ern_genes)
        #    print("[-ADDING ERN GENES TO SOLVE-RD EXOMEDEPTH CALLS-]")
        #    add_ern_genes_to_cnvcalls(exd_calls, ern_genes)
        #    compare_gatk_exomedepth_erngenes()


    # Perform the comparison with the Conifer calls
    if compare_params["callstool"] == "conifer":
        print("[-READING SOLVE-RD CONIFER CALLS-]")
        conifer_calls = read_conifer(compare_params["calls"])

        print("[-COMPARING GATK4 AND SOLVE-RD CONIFER CALLS-]")
        compare_gatk_solverd_overlap(ccrs_calls, conifer_calls, compare_params["percentage-overlap"], CONIFER_TO_GATK4)
        unique_gatk_samples = set(ccrs_calls.keys()) - set(conifer_calls.keys())
        unique_conifer_samples = set(conifer_calls.keys()) - set(ccrs_calls.keys())

        print("[-WRITING CONIFER CALLS WITH WHICH CCRS CALLS OVERLAP TO FILE-]")
        wrote_file = write_overlapping_calls(f"{outpath}_overlapping.txt", conifer_calls, "Conifer")
        print(f"...Wrote overlapping calls to file?: {wrote_file}...")

        print("[-WRITING UNIQUE GATK4 AND CONIFER SAMPLE NAMES TO FILE-]")
        wrote_gatk4 = write_unique_samples(f"{outpath}_unique_gatk_samples.txt", unique_gatk_samples)
        print(f"...Wrote unique GATK4 sample names?: {wrote_gatk4}...")
        wrote_conifer = write_unique_samples(f"{outpath}_unique_conifer_samples.txt", unique_conifer_samples)
        print(f"...Wrote unique Conifer sample names?: {wrote_conifer}...")

        # Display the numbers
        count_samples_calls_without_ccrs(conifer_calls)


    # Perform the comparison with the ClinCNV calls
    if compare_params["callstool"] == "clincnv":
        print("[-READING SOLVE-RD CLINCNV CALLS-]")
        clincnv_calls = read_clincnv(compare_params["calls"])

        print("[-COMPARING GATK4 AND SOLVE-RD CLINCNV CALLS-]")
        compare_gatk_solverd_overlap(ccrs_calls, clincnv_calls, compare_params["percentage-overlap"], CLINCNV_TO_GATK4)
        unique_gatk_samples = set(ccrs_calls.keys()) - set(clincnv_calls.keys())
        unique_clincnv_samples = set(clincnv_calls.keys()) - set(ccrs_calls.keys())

        print("[-WRITING CLINCNV CALLS WITH WHICH CCRS CALLS OVERLAP TO FILE-]")
        wrote_file = write_overlapping_calls(f"{outpath}_overlapping.txt", clincnv_calls, "ClinCNV")
        print(f"...Wrote overlapping calls to file?: {wrote_file}...")

        print("[-WRITING UNIQUE GATK4 AND CLINCNV SAMPLE NAMES TO FILE-]")
        wrote_gatk4 = write_unique_samples(f"{outpath}_unique_gatk_samples.txt", unique_gatk_samples)
        print(f"...Wrote unique GATK4 sample names?: {wrote_gatk4}...")
        wrote_clincnv = write_unique_samples(f"{outpath}_unique_clincnv_samples.txt", unique_clincnv_samples)
        print(f"...Wrote unique ClinCNV sample names?: {wrote_clincnv}...")

        # Dispaly the numbers
        count_samples_calls_without_ccrs(clincnv_calls)


    # Perform the comparison with the VarGenius calls
    if compare_params["callstool"] == "vargenius":
        print("[-READING SOLVE-RD VARGENIUS CALLS-]")
        vargenius_calls = read_vargenius(compare_params["calls"])

        print("[-COMPARING GATK4 CALLS AND SOLVE-RD VARGENIUS CALLS-]")
        compare_gatk_solverd_overlap(ccrs_calls, vargenius_calls, compare_params["percentage-overlap"], VARGENIUS_TO_GATK)
        unique_gatk_samples = set(ccrs_calls.keys()) - set(vargenius_calls.keys())
        unique_vargenius_samples = set(vargenius_calls.keys()) - set(ccrs_calls.keys())

        print("[-WRITING VARGENIUS CALLS WITH WHICH CCRS CALLS OVERLAP TO FILE-]")
        wrote_file = write_overlapping_calls(f"{outpath}_overlapping.txt", vargenius_calls, "VarGenius")
        print(f"...Wrote overlapping calls to file?: {wrote_file}...")

        print("[-WRITING UNIQUE GATK4 AND CLINCNV SAMPLE NAMES TO FILE-]")
        wrote_gatk4 = write_unique_samples(f"{outpath}_unique_gatk_samples.txt", unique_gatk_samples)
        print(f"...Wrote unique GATK4 sample names?: {wrote_gatk4}...")
        wrote_vargenius = write_unique_samples(f"{outpath}_unique_vargenius_samples.txt", unique_vargenius_samples)
        print(f"...Wrote unique VarGenius sample names?: {wrote_vargenius}...")

        # Display the numbers
        count_samples_calls_without_ccrs(vargenius_calls)


def tmp_show_exd_overlaps(exdcalls):
    """."""
    good_overlap = 0
    bad_overlap = 0
    for samplename in exdcalls:
        for chromname in exdcalls[samplename]:
            for exdcall in exdcalls[samplename][chromname]:
                # print(f"Overlapping calls for EXD {exdcall.sample_name}_{exdcall.cnv_chrom}:{exdcall.cnv_start}-{exdcall.cnv_end} ({exdcall.cnv_call})")
                if len(exdcall.overlapping_ccrs) > 0:
                    print(f"Overlapping calls for EXD {exdcall.sample_name}_{exdcall.cnv_chrom}:{exdcall.cnv_start}-{exdcall.cnv_end} ({exdcall.cnv_call})")
                    good_overlap += 1
                    for ccrscall in exdcall.overlapping_ccrs:
                        print(f"\t{ccrscall.ccrs_sample}_{ccrscall.get_region_string()} ({ccrscall.ccrs_call})")
                if len(exdcall.overlapping_ccrs_2) > 0:
                    print(f"Overlapping calls for EXD {exdcall.sample_name}_{exdcall.cnv_chrom}:{exdcall.cnv_start}-{exdcall.cnv_end} ({exdcall.cnv_call})")
                    bad_overlap += 1
                    for ccrscall in exdcall.overlapping_ccrs_2:
                        print(f"\t{ccrscall.ccrs_sample}_{ccrscall.get_region_string()} ({ccrscall.ccrs_call})")
    print(f"Good overlap: {good_overlap}")
    print(f"Bad overlap: {bad_overlap}")


def tmp_show_unique_samples(ugatksamples, uothersamples, otherlabel):
    """Show unique samples for gatk4 and the other tool."""
    print(f"...Unique GATK4 samples ({len(ugatksamples)})...")
    print(ugatksamples)
    print(f"...Unique {otherlabel} samples ({len(uothersamples)})...")
    print(uothersamples)


def write_overlapping_calls(outfileloc, solverdcallswithccrs, labelname):
    """Write CCRS calls that overlap with other Solve-RD calls.

    Parameters
    ----------
    outfileloc : str
        Path to write output file to
    solverdcallswithccrs : dict
        Solve-RD calls with which CCRS calls overlap >=80%
    """
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(f"Sample\t{labelname}_CNV\t{labelname}_Call\tCCRS_CNV\tCCRS_Call\n")
            for samplename in solverdcallswithccrs:
                for chromname in solverdcallswithccrs[samplename]:
                    for exdcall in solverdcallswithccrs[samplename][chromname]:
                        if len(exdcall.overlapping_ccrs) > 0:
                            for ccrscall in exdcall.overlapping_ccrs:
                                outfile.write(f"{samplename}\t{exdcall.get_region_string()}\t{exdcall.cnv_call}\t{ccrscall.get_region_string()}\t{ccrscall.ccrs_call}\n")
                        if len(exdcall.overlapping_ccrs_2) > 0:
                            for ccrscall in exdcall.overlapping_ccrs_2:
                                outfile.write(f"{samplename}\t{exdcall.get_region_string()}\t{exdcall.cnv_call}\t{ccrscall.get_region_string()}\t{ccrscall.ccrs_call}\n")
        file_written = True
    except IOError:
        print("Could not write calls to file")
    finally:
        return file_written


def write_unique_samples(outfileloc, samplenames):
    """Write unique samples to file.

    Parameters
    ----------
    outfileloc : str
        Path to write output file
    samplenames : set
        Unique sample names
    """
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\n")
            for samplename in samplenames:
                outfile.write(f"{samplename}\n")
        file_written = True
    except IOError:
        print("Could not write unique sample names to file")
    finally:
        return file_written


def show_overlapping_calls(outfileloc, solverdcallswithccrs, labelname):
    print(f"Sample\t{labelname}_CNV\t{labelname}_Call\tCCRS_CNV\tCCRS_Call\n")
    for samplename in solverdcallswithccrs:
        for chromname in solverdcallswithccrs[samplename]:
            for exdcall in solverdcallswithccrs[samplename][chromname]:
                if len(exdcall.overlapping_ccrs) > 0:
                    for ccrscall in exdcall.overlapping_ccrs:
                        print(f"{samplename}\t{exdcall.get_region_string()}\t{exdcall.cnv_call}\t{ccrscall.get_region_string()}\t{ccrscall.ccrs_call}\n")
                if len(exdcall.overlapping_ccrs_2) > 0:
                    for ccrscall in exdcall.overlapping_ccrs_2:
                        print(f"{samplename}\t{exdcall.get_region_string()}\t{exdcall.cnv_call}\t{ccrscall.get_region_string()}\t{ccrscall.ccrs_call}\n")


def count_samples_calls_without_ccrs(solverdcalls):
    """Show samples and ."""
    totalcalls = 0
    totalsamples = 0
    nosamplecount = 0
    nocallcount = 0
    hassamplecount = 0
    hascallcount = 0

    for samplename in solverdcalls:
        totalsamples += 1
        hasccrscount = 0
        for chromname in solverdcalls[samplename]:
            for solverdcall in solverdcalls[samplename][chromname]:
                totalcalls += 1
                if len(solverdcall.overlapping_ccrs) == 0 and len(solverdcall.overlapping_ccrs_2) == 0:
                    callcount += 1
                else:
                    hasccrscount += 1
                    hascallcount += 1
        if hasccrscount == 0:
            nosamplecount += 1
        else:
            hassamplecount += 1

    # Display the counts
    print(f"Total Solve-RD samples: {totalsamples}")
    print(f"Total Solve-RD calls: {totalcalls}")
    print(f"Solve-RD samples with overlapping GATK4 calls: {hassamplecount}/{totalsamples} ({round((hassamplecount/totalsamples)*100, 2)}%)")
    print(f"Solve-RD calls with overlapping GATK4 calls: {hascallcount}/{totalcalls} ({round((hascallcount/totalcalls)*100, 2)}%)")
    print(f"Solve-RD samples without overlapping GATK4 calls: {nosamplecount}/{totalsamples} ({round((nosamplecount/totalsamples)*100, 2)}%)")
    print(f"Solve-RD calls without overlapping GATK4 calls: {nocallcount}/{totalcalls} ({round((nocallcount/totalcalls)*100, 2)}%)")


if __name__ == "__main__":
    main()
