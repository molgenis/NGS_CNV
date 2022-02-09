import argparse


def get_params():
    fil_args = argparse.ArgumentParser()
    fil_args.add_argument("-l", "--intervallist", required=True, dest="intervallist", help="Path to intervallist to filter")
    fil_args.add_argument("-c", "--chromosome", required=True, dest="chromosome", help="Chromosome to filter on")
    fil_args.add_argument("-o", "--outfile", required=True, dest="outfile", help="Path to filtered intervallist to")
    return vars(fil_args.parse_args())


def main():
    fil_params = get_params()
    try:
        outfile = open(fil_params["outfile"], 'w')
        with open(fil_params["intervallist"], 'r') as intervallist:
            for fileline in intervallist:
                if fileline.startswith("@"):
                    outfile.write(fileline)
                else:
                    if fileline.startswith(fil_params["chromosome"]):
                        outfile.write(fileline)
        outfile.close()
    except IOError:
        print("Could not filter interval list file")


if __name__ == "__main__":
    main()
