#!/usr/bin/env python
import argparse
import gzip


def get_params():
    cacl_args = argparse.ArgumentParser()
    cacl_args.add_argument("-g", "--gnomad", required=True, dest="gnomad", help="Path to gnomAD VCF file")
    cacl_args.add_argument("-o", "--outfile", required=True, dest="outfile", help="Path to write filtered output file to")
    cacl_args.add_argument("-f", "--allele-frequency", dest="allele-frequency", default=10, type=int, help="Allele Frequency to filter on")
    cacl_args.add_argument("-a", "--allele-field", dest="allele-field", default="AF", help="Allele info field to use")
    return vars(cacl_args.parse_args())


def get_allele_frequency(gnomadinfo, affield):
    for infofield in gnomadinfo:
        infodata = infofield.split("=")
        if infodata[0] == affield:
            return infodata[1]


def filter_gnomad(gnomadloc, outfileloc, af_cutoff, af_field):
    try:
        outfile = open(outfileloc, 'w')

        # Start filtering the gnomAD data
        with gzip.open(gnomadloc, 'rt') as gnomadfile:
            for gnomadline in gnomadfile:
                if gnomadline.startswith('#'):
                    outfile.write(gnomadline)
                else:
                    gnomaddata = gnomadline.strip().split()
                    gnomadinfo = gnomaddata[7].split(";")
                    allele_frequency = get_allele_frequency(gnomadinfo, af_field)

                    # Check that the allele frequency is equal or larger than the set cutoff
                    if allele_frequency != ".":
                        if float(allele_frequency)*100 >= af_cutoff:
                            outfile.write(gnomadline)
        outfile.close()
    except IOError:
        print("")

def main():
    cacl_params = get_params()
    filter_gnomad(cacl_params["gnomad"], cacl_params["outfile"], cacl_params["allele-frequency"], cacl_params["allele-field"])


if __name__ == "__main__":
    main()
