import sys

def read_ccrs_calls(ccrsfileloc):
    """Read the CCRS data and return the numbers of probes per sample."""
    ccrs_calls = {}
    try:
        with open(ccrsfileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[0] not in ccrs_calls:
                    ccrs_calls[filelinedata[0]] = []
                ccrs_calls[filelinedata[0]].append(f"{filelinedata[1]}:{filelinedata[2]}-{filelinedata[3]}")
    except IOError:
        print("Could not read combined CCRS file")
    finally:
        return ccrs_calls


def write_num_of_calls(outfileloc, ccrscalls):
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("Sample\tNum_Calls\n")
            for samplename in ccrscalls:
                outfile.write(f"{samplename}\t{len(ccrscalls[samplename])}\n")
        file_written = True
    except IOError:
        print("Could not write num of calls file")
    finally:
        return file_written


def main():
    ccrs_data = read_ccrs_calls(sys.argv[1])
    wrote_file = write_num_of_calls(sys.argv[2], ccrs_data)
    print(f"Wrote output file?: {wrote_file}")


if __name__ == "__main__":
    main()
