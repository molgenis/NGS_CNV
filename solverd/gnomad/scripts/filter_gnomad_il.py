import argparse
import gzip


def get_params():
    cacl_args = argparse.ArgumentParser()
    cacl_args.add_argument("-g", "--gnomad", required=True, dest="gnomad", help="Path to gnomAD VCF file")
    cacl_args.add_argument("-l", "--intervallist", required=True, dest="intervallist", help="Path to preprocess intervallist")
    cacl_args.add_argument("-o", "--outfile", required=True, dest="outfile", help="Path to write filtered output file to")
    cacl_args.add_argument("-H", "--header", dest="header", action="store_true", help="Write the header?")
    return vars(cacl_args.parse_args())


def read_interval_list(intervalfileloc):
    interval_data = {}
    try:
        with open(intervalfileloc, 'r') as intervallist:
            for fileline in intervallist:
                if not fileline.startswith("@"):
                    intervaldata = fileline.strip().split()
                    if intervaldata[0] not in interval_data:
                        interval_data[intervaldata[0]] = []
                    interval_data[intervaldata[0]].append(range(int(intervaldata[1]), int(intervaldata[2])+1))
    except IOError:
        print("")
    finally:
        return interval_data


def gnomad_snp_in_interval(chromintervals, gnomadpos):
    for chromintv in chromintervals:
        if gnomadpos in chromintv:
            return True
    return False


def filter_gnomad(gnomadfileloc, interval_data, write_header, outfileloc):
    try:
        outfile = open(outfileloc, 'w')

        # Start filtering the gnomad VCF file
        with gzip.open(gnomadfileloc, 'rt') as gnomadfile:
            for gnomadline in gnomadfile:
                if gnomadline.startswith('#') and write_header:
                    outfile.write(gnomadline)
                else:
                    gnomaddata = gnomadline.strip().split()
                    if gnomaddata[0] in interval_data:
                        if gnomad_snp_in_interval(interval_data[gnomaddata[0]], gnomaddata[1]):
                            outfile.write(gnomadline)
        outfile.close()
    except IOError:
        print("")


def main():
    cacl_params = get_params()
    intervallistdata = read_interval_list(cacl_params["intervallist"])
    filter_gnomad(cacl_params["gnomad"], intervallistdata, cacl_params["header"], cacl_params["outfile"])


if __name__ == "__main__":
    main()
