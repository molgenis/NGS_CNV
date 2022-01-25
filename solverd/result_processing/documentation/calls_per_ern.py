import sys

def read_rd3_samples_to_ern(infileloc):
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
    # ernlabel = "ERN EURO-NMD"
    # ernlabel = "ERN GENTURIS"
    # ernlabel = "ERN ITHACA"
    # ernlabel = "ERN-RND"
    # ernlabel = "UDN-Spain"
    ernlabel = sys.argv[4]
    ern_samples = read_rd3_samples_to_ern(sys.argv[1])
    ern_counts = get_calls_for_ern(sys.argv[2])
    zero_count_samples = set(ern_samples[ernlabel]) - set(ern_counts.keys())
    sample_names = ern_samples[ernlabel]
    sample_names.sort()
    wrote_file = write_sample_counts(sys.argv[3], sample_names, zero_count_samples, ern_counts)
    print(f"...{wrote_file}...")


if __name__ == "__main__":
    main()
