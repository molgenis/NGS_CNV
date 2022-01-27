# import sys
import argparse


def get_params():
    """Define, receive and return set CLI parameter values."""
    ERN_CHOICES = ["ERN EURO-NMD", "ERN GENTURIS", "ERN ITHACA", "ERN-RND", "UDN-Spain"]
    cpe_args = argparse.ArgumentParser()
    cpe_args.add_argument("-i", "--infile", type=str, required=True, dest="infile", help="Path to GATK4 CNV calls per ERN")
    cpe_args.add_argument("-s", "--samples-to-ern", type=str, required=True, dest="samples-to-ern", help="Path to samples to ERN file")
    cpe_args.add_argument("-o", "--outfile", type=str, required=True, dest="outfile", help="Path to write output file to")
    cpe_args.add_argument("-l", "--label", type=str, required=True, dest="label", choices=ERN_CHOICES, help="Name of the ERN")
    return vars(cpe_args.parse_args())


def read_rd3_samples_to_ern(infileloc):
    """Read the samples to ERN file.

    Parameters
    ----------
    infileloc : str
        Path to the file with samples to ERN.

    Returns
    -------
    ernsamples : dict
        Sample names per ERN
    """
    ernsamples = {}
    try:
        with open(infileloc, 'r') as infile:
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[1] not in ernsamples:
                    ernsamples[filelinedata[1]] = []
                ernsamples[filelinedata[1]].append(filelinedata[0])
    except IOError:
        print("aap")
    finally:
        return ernsamples


def get_calls_for_ern(infileloc):
    """Count the occurrence (calls) of each sample.

    Parameters
    ----------
    infileloc : str
        Path to GATK4 calls per ERN

    Returns
    -------
    erncallcounts : dict
        Call counts per sample name
    """
    erncallcounts = {}
    try:
        with open(infileloc, 'r') as infile:
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[0] not in erncallcounts:
                    erncallcounts[filelinedata[0]] = 0
                erncallcounts[filelinedata[0]] += 1
    except IOError:
        print("noot")
    finally:
        return erncallcounts


def write_sample_counts(outfileloc, samplenames, zerocountsamples, countsamples):
    """Write the sample call counts to a tab separated file.

    Parameters
    ----------
    outfileloc : str
        Path to write output file to
    samplenames : list of str
        Sorted list of sample names
    zerocountsamples : list of str
        List of sample names without calls
    countsamples : dict
        Samples and their call counts
    """
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tCount\n")
            for samplename in samplenames:
                if samplename in zerocountsamples:
                    outfile.write(f"{samplename}\t0\n")
                else:
                    outfile.write(f"{samplename}\t{countsamples[samplename]}\n")
        file_written = True
    except IOError:
        print("mies")
    finally:
        return file_written


def main():
    cpe_params = get_params()
    ernlabel = cpe_params["label"]
    ern_samples = read_rd3_samples_to_ern(cpe_params["samples-to-ern"])
    ern_counts = get_calls_for_ern(cpe_params["infile"])
    zero_count_samples = set(ern_samples[ernlabel]) - set(ern_counts.keys())
    sample_names = ern_samples[ernlabel]
    sample_names.sort()
    wrote_file = write_sample_counts(cpe_params["outfile"], sample_names, zero_count_samples, ern_counts)
    print(f"...Wrote output file?: {wrote_file}...")


if __name__ == "__main__":
    main()
