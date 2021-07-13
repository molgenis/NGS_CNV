import argparse
import gzip


def get_params():
    gnomad_args = argparse.ArgumentParser()
    gnomad_args.add_argument("-g", "--gnomad", required=True, dest="gnomad", help="")
    gnomad_args.add_argument("-o", "--outfile", required=True, dest="outfile", help="")
    return vars(gnomad_args.parse_args())


def main():
    gnomad_params = get_params()
    try:
        outfile = open(gnomad_params["outfile"], 'w')
        with gzip.open(gnomad_params["gnomad"], "rt") as gnomadfile:
            for gnomadline in gnomadfile:
                if gnomadline.startswith("#"):
                    outfile.write(gnomadline)
                else:
                    gnomaddata = gnomadline.strip().split()
                    varref = gnomaddata[3].split(",")
                    varalt = gnomaddata[4].split(",")
                    if len(max(varref, key=len)) == 1 and len(max(varalt, key=len)) == 1:
                        outfile.write(gnomadline)
        outfile.close()
    except IOError:
        print("Could not process gnomAD file")


if __name__ == "__main__":
    main()
