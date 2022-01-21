import os
import argparse

def get_params():
    """Define, receive and return set parameter values."""
    fn_args = argparse.ArgumentParser()
    fn_args.add_argument("-i", "--indir", type=str, dest="indir", required=True, help="Path to CCRS directory to filter neutrals from")
    fn_args.add_argument("-o", "--outdir", type=str, dest="outdir", required=True, help="Path to output directory")
    fn_args.add_argument("-de", "--deletion-outliers", type=str, dest="deletion-outliers", help="Path to file with deletion outliers")
    fn_args.add_argument("-du", "--duplication-outliers", type=str, dest="duplication-outliers", help="Path to file with duplication outliers")
    return vars(fn_args.parse_args())


def read_outliers(outlierfileloc):
    """Reads a file containing a list of outlier samples."""
    outliersamples = []
    try:
        with open(outlierfileloc, 'r') as outlierfile:
            for fileline in outlierfile:
                outliersamples.append(fileline.strip())
    except IOError:
        print("Could not read outliers file")
    finally:
        return outliersamples


def filter_neutrals_from_file(infileloc, outfileloc):
    """Filter neutrals from the CCRS input file and write the filtered data to a new output file."""
    filter_success = False
    try:
        outfile = open(outfileloc, 'w')

        # Start filtering the input file and write the results to a new output file to retain the original
        with open(infileloc, 'r') as infile:
            outfile.write(next(infile))
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[5] != '0':
                    outfile.write(fileline)

        outfile.close()
        filter_success = True
    except IOError:
        print("")
    finally:
        return filter_success


def main():
    """Do the main work."""
    fn_params = get_params()

    # Form the input and output directory paths
    indir = fn_params["indir"]
    indir = f"{indir}/" if not indir.endswith("/") else indir
    outdir = fn_params["outdir"]
    outdir = f"{outdir}/" if not outdir.endswith("/") else outdir

    # Gather the CCRS files
    ccrs_files = [x for x in os.listdir(indir) if x.endswith(".igv.seg")]

    # Set the lists of outliers (will be empty if no file is provided)
    deletion_outliers = read_outliers(fn_params["deletion-outliers"]) if fn_params["deletion-outliers"] else []
    duplication_outliers = read_outliers(fn_params["duplication-outliers"]) if fn_params["duplication-outliers"] else []

    # Start filtering neutrals from the CCRS files
    unfiltered_files = []
    for ccrsfile in ccrs_files:
        file_was_filtered = False
        samplename = ccrsfile.split(".")[0]
        if samplename not in deletion_outliers and samplename not in duplication_outliers:
            file_was_filtered = filter_neutrals_from_file(f"{indir}{ccrsfile}", f"{outdir}{ccrsfile}")
        else:
            print(f"...Sample {samplename} is an outlier and will therefore be skipped...")

        if not file_was_filtered:
            unfiltered_files.append(f"{indir}{ccrsfile}")

    # Display the unfiltered files, if there are any
    if len(unfiltered_files) > 0:
        print(f"Unfiltered samples: {unfiltered_files}")


if __name__ == "__main__":
    main()
