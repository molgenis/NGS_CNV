import argparse
from ccrscall import CcrsCall
from gnomadentry import GnomadEntry


def get_params():
    """Define, receive and return set CLI parameter values."""
    gnomad_filter_args = argparse.ArgumentParser()
    gnomad_filter_args.add_argument("-a", "--annotate", action="store_true", dest="annotate", help="Only annotate results")
    gnomad_filter_args.add_argument("-af", "--annotation-file", type=str, dest="annotation-file", help="Path to existing gnomAD annotation file")
    gnomad_filter_args.add_argument("-f", "--frequency", type=float, dest="frequency", help="Miinmum frequency for filtering")
    gnomad_filter_args.add_argument("-ff", "--filter-field", type=str, dest="filter-field", help="")
    gnomad_filter_args.add_argument("-g", "--gnomad-file", type=str, dest="gnomad-file", help="Path to gnomAD file")
    gnomad_filter_args.add_argument("-i", "--infile", type=str, dest="infile", required=True, help="Path to combined CCRS input file to")
    gnomad_filter_args.add_argument("-o", "--outdir", type=str, dest="outdir", required=True, help="Path to output directory")
    gnomad_filter_args.add_arugment("-p", "--output-prefix", type=str, dest="output-prefix", help="Output prefix to use for files")
    return vars(gnomad_filter_args.parse_args())


def read_combined_ccrs(infileloc):
    """Read the combined CCRS file."""
    ccrs_calls = {}
    try:
        with open(infileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")

                # Add the sample to the CCRS calls data
                if filelinedata[0] not in ccrs_calls:
                    ccrs_calls[filelinedata[0]] = []
                ccrs_calls[filelinedata[0]].append(CcrsCall(filelinedata[0], filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), int(filelinedata[4]), filelinedata[5], float(filelinedata[6])))
    except IOError:
        print("Could not read combined CCRS file")
    finally:
        return ccrs_calls


def read_annotation_file(infileloc):
    """."""
    annotationdata = {}
    try:
        with open(infileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
    except IOError:
        print("Could not write annotation file")
    finally:
        return annotationdata


def read_gnomad_file(infileloc):
    """Read the gnomAD file."""
    gnomaddata = {}
    try:
        with open(infileloc, 'r') as infile:
            
    except IOError:
        print("Could not read ")
    finally:
        return gnomaddata


def write_annotation_file(annotatedccrsdata, samplenames, outfileloc):
    """."""
    wrote_file = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Samplet\tChromosome\tStart\tEnd\tNum_Probes\tCall\tSegment_Mean\t\n")
            for samplename in samplenames:
                for accrscall in annotatedccrsdata[samplename]:
                    outfile.write(f"{accrscall.ccrs_chrom}\t{accrscall.ccrs_start}\t{accrscall.ccrs_end}\t{accrscall.ccrs_numofprobes}\t{}\t{}\t{}")
        wrote_file = True
    except IOError:
        print("")
    finally:
        return wrote_file


def annotate_ccrs(ccrsdata, gnomadentries):
    """Annotate the CCRS data with gnomAD entries."""


def main():
    """Either annotate or filter CCRS calls with gnomAD."""
    gnomad_filter_params = get_params()
    outdir = gnomad_filter_params["outdir"]+"/" if not gnomad_filter_params["outdir"].endswith("/") else gnomad_filter_params["outdir"]
    outprefix = gnomad_filter_params["output-prefix"]
    ccrs_data = read_combined_ccrs(gnomad_filter_params["infile"])

    sample_names = list(ccrs_data.keys())
    sample_names.sort()

    # Determine whether to annotate or filter
    if gnomad_filter_params["annotate"]:
        #Do annotation
        print("[-START ANNOTATING CCRS CALLS WITH GNOMAD DATA-]")
        if gnomad_filter_params["annotation-file"]:
            print("...gnomAD annotation file has been provided, so the CCRS calls have already been annotated...")
        else:
            annotate_ccrs()
            write_annotation_file(f"{outdir}{outprefix}.annotated.txt")
    else:
        #Do filtering
        if gnomad_filter_params["annotation-file"]:
            read_annotation_file()


if __name__ == "__main__":
    main()
