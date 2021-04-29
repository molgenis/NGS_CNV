import argparse


def get_params():
    gnomad_args = argparse.ArgumentParser()
    gnomad_args.add_argument("-i", "--infile", required=True, dest="infile", help="Path to interval list to convert")
    gnomad_args.add_argument("-o", "--outfile", required=True, dest="outfile", help="Path to write BED output file to")
    gnomad_args.add_argument("-H", "--header", dest="header", action="store_true", default=True, help="Give BED file a header?")
    gnomad_args.add_argument("-c", "--chromosomes", dest="chromosomes", nargs="+", default="*", help="")
    return vars(gnomad_args.parse_args())


def main():
    gnomad_params = get_params()
    include_all_chroms = "*" in gnomad_params["chromosomes"]
    try:
        bedfile = open(gnomad_params["outfile"], 'w')
        if gnomad_params["header"]:
            bedfile.write("chrom\tstart\tend\n")

        with open(gnomad_params["infile"], 'r') as intervallist:
            for intervalline in intervallist:
                if not intervalline.startswith("@"):
                    if include_all_chroms:
                        bedfile.write(intervalline)
                    else:
                        intervaldata = intervalline.strip().split()
                        if intervaldata[0] in gnomad_params["chromosomes"]:
                            bedfile.write(intervalline)
        bedfile.close()


if __name__ == "__main__":
    main()
