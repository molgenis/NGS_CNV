import argparse



def get_params():
    """Define, receive and return parameter values"""
    dupdel_args = argparse.ArgumentParser()
    dupdel_args.add_argument("-i", "--infile", type=str, required=True, dest="infile", help="")
    dupdel_args.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", help="")
    dupdel_args.add_argument("-p", "--prefix", type=str, required=True, dest="prefix", help="")
    dupdel_args.add_argument("-m1", "--min-duplication-outlier", type=int, required=True, dest="min-duplication-outlier")
    dupdel_args.add_argument("-m2", "--min-deletion-outlier", type=int, required=True, dest="min-deletion-outlier")
    return vars(dupdel_args.parse_args())


def filter_outliers(infileloc, outfileloc, mindup, mindel):
    """Filter outlier samples based on number of duplications and deletions"""
    dup_outlier_samples = []
    del_outlier_samples = []
    try:
        outfile = open(outfileloc, 'w')

        with open(infileloc, 'r') as infile:
            outfile.write(next(infile))
            for fileline in infile:
                dupdeldata = fileline.strip().split("\t")
                write_line = True

                if int(dupdeldata[1]) >= mindup:
                    dup_outlier_samples.append(dupdeldata[0])
                    write_line = False
                if int(dupdeldata[2]) >= mindel:
                    del_outlier_samples.append(dupdeldata[0])
                    write_line = False

                if write_line:
                    outfile.write(fileline)
        outfile.close()
    except IOError:
        print("Could not filter outliers")
    finally:
        return {"dup":dup_outlier_samples, "del":del_outlier_samples}


def write_outliers(outliersfileloc, outliersdata):
    """Write the filtered duplication or deletion outliers to file"""
    try:
        with open(outliersfileloc, 'w') as outliersfile:
            # Write the outlier samples to file
            for outliersample in outliersdata:
                outliersfile.write(f"{outliersample}\n")
            outliersfile.write("\n")
    except IOError:
        print("Could not write outlier samples to file")


def main():
    """Do the main work"""
    dupdel_params = get_params()

    outdir = dupdel_params["outdir"]
    outdir = f"{outdir}/" if not outdir.endswith("/") else outdir
    outprefix = dupdel_params["prefix"]

    batch_outliers = filter_outliers(dupdel_params["infile"], f"{outdir}{outprefix}_filtered_dupdel.txt", dupdel_params["min-duplication-outlier"], dupdel_params["min-deletion-outlier"])
    write_outliers(f"{outdir}{outprefix}_dup_outliers.txt", batch_outliers["dup"])
    write_outliers(f"{outdir}{outprefix}_del_outliers.txt", batch_outliers["del"])


if __name__ == "__main__":
    main()
