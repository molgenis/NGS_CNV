import os
import argparse

def get_params():
    """Define, receive and return parameters values."""
    combccrs_args = argparse.ArgumentParser()
    combccrs_args.add_argument("-i", "--indir", type=str, dest="indir", required=True, help="Path ")
    combccrs_args.add_argument("-o", "--outfile", type=str, dest="outfile", required=True, help="Path to write output file to")
    return vars(combccrs_args.parse_args())


def combine_ccrs(outfileloc, indirloc, ccrsfiles, samplenames):
    """Combines the CCRS file into one."""
    wrote_file = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tChromosome\tStart\tEnd\tNum_Probes\tCall\tSegment_Mean\n")
            for samplename in samplenames:
                print(f"...Adding {ccrsfiles[samplename]}...")
                with open(f"{indirloc}{ccrsfiles[samplename]}", 'r') as ccrsfile:
                    next(ccrsfile)
                    for fileline in ccrsfile:
                        outfile.write(fileline)
        wrote_file = True
    except IOError:
        print("Something went wrong combining the CCRS fiels into one")
    finally:
        return wrote_file


def main():
    """Do the main work."""
    combccrs_params = get_params()

    # Make sure the path to the indir ends with a /
    indir = combccrs_params["indir"]
    indir = f"{indir}/" if not indir.endswith("/") else indir

    # Collect the input files and sample names
    ccrs_files = {x.split(".")[0]:x for x in os.listdir(indir) if x.endswith(".igv.seg")}
    sample_names = list(ccrs_files.keys())
    sample_names.sort()

    # Combine the CCRS files
    wrote_combinedfile = combine_ccrs(combccrs_params["outfile"], indir, ccrs_files, sample_names)
    print(f"Wrote combined CCRS file?: {wrote_combinedfile}")


if __name__ == "__main__":
    main()
