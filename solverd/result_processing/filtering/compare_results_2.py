#!/usr/bin/env python
import argparse
from ccrscall import CcrsCall
from read_combined_ccrs import read_combined_ccrs
from solverdcall import SolveRdCall
from solverdexomedepthcall import SolveRdExomeDepthCall


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
    """."""
    compare_args = argparse.ArgumentParser()
    compare_args.add_argument("-c", "--clincnv", type=str, required=True, dest="clincnv", help="")
    compare_args.add_argument("-f", "--conifer", type=str, required=True, dest="conifer", help="")
    compare_args.add_argument("-e", "--exomedepth", type=str, required=True, dest="exomedepth", help="")
    compare_args.add_argument("-g", "--gatk4", type=str, required=True, dest="gatk4", help="")
    compare_args.add_argument("-v", "--vargenius", type=str, required=True, dest="vargenius", help="")
    compare_args.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", help="")
    compare_args.add_argument("-op", "--out-prefix", type=str, required=True, dest="out-prefix", help="")
    return vars(compare_args.parse_args())


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

                if filelinedata[7] not in conifer_calls:
                    conifer_calls[filelinedata[7]] = {}
                if filelinedata[1] not in conifer_calls[filelinedata[7]]:
                    conifer_calls[filelinedata[7]][filelinedata[1]] = []
                conifer_calls[filelinedata[7]][filelinedata[1]].append(SolveRdCall(filelinedata[7], filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), filelinedata[4], filelinedata[5]))
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

                if filelinedata[0] not in vargenius_calls:
                    vargenius_calls[filelinedata[0]] = {}
                if filelinedata[13] not in vargenius_calls[filelinedata[0]]:
                    vargenius_calls[filelinedata[0]][filelinedata[13]] = []
                vargenius_calls[filelinedata[0]][filelinedata[13]].append(SolveRdCall(filelinedata[0], filelinedata[13], int(filelinedata[14]), int(filelinedata[15]), int(filelinedata[18]), filelinedata[16]))
    except IOError:
        print("Could not read VarGenius file")
    finally:
        return vargenius_calls


def main():
    """."""
    compare_params = get_params()
    compare_params["outdir"] = compare_params["outdir"]+"/" if not compare_params["outdir"].endswith("/") else compare_params["outdir"]
    outpath = compare_params["outdir"] + compare_params["out-prefix"]

    # Read the Solve-RD CNV calls files
    print("[-READING SOLVERD CLINCNV CNV CALLS-]")
    clincnv_calls = read_clincnv(compare_params["clincnv"])
    print("[-READING SOLVERD CONIFER CNV CALLS-]")
    conifer_calls = read_conifer(compare_params["conifer"])
    print("[-READING SOLVERD EXOMEDEPTH CNV CALLS-]")
    exomedepth_calls = read_exomedepth(compare_params["exomedepth"])
    print("[-READING SOLVERD VARGENIUS CNV CALLS-]")
    vargenius_calls = read_vargenius(compare_params["vargenius"])

    # Read the GATK4 CNV calls
    print("[-READING GATK4 CNV CALLS-]")
    ccrs_data = read_combined_ccrs(compare_params["gatk4"])
    ccrs_header = ccrs_data[0]
    gatk4_calls = ccrs_data[1]

    # print(clincnv_calls.keys())
    # print(conifer_calls.keys())
    # print(exomedepth_calls.keys())
    # print(gatk4_calls.keys())
    # print(vargenius_calls.keys())

    # Determine the shared samples between tools
    print("[-DETERMENING THE SAMPLES SHARED BETWEEN ALL TOOLS-]")
    shared_samples_all = set(clincnv_calls.keys()) & set(conifer_calls.keys()) & set(exomedepth_calls.keys()) & set(vargenius_calls.keys()) & set(gatk4_calls.keys())
    print(f"Samples shared between all tools: {len(shared_samples_all)}")


    # shared_samples_gatk4clincnv = set(gatk4_calls.keys()) & set(clincnv_calls.keys())
    # shared_samples_gatk4conifer = set(gatk4_calls.keys()) & set(conifer_calls.keys())
    # shared_samples_gatk4exomedepth = set(gatk4_calls.keys()) & set(exomedepth_calls.keys())
    # shared_samples_gatk4vargenius = set(gatk4_calls.keys()) & set(vargenius_calls.keys())
    # shared_samples_ = set(.keys()) & set(.keys())
    # shared_samples_ = set(.keys()) & set(.keys())


if __name__ == "__main__":
    main()


# Usage:
# python compare_results_2.py
# -c /path/to/clincnv_calls.tsv
# -e /path/to/exomedepth_calls.tsv
# -f /path/to/conifer_calls.tsv
# -g /path/to/gatk4_calls.tsv
# -v /path/to/vargenius_calls.tsv
# -o /path/to/outdir/
# -op compare_all
