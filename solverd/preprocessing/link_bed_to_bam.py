import argparse


def get_params():
    lb2b_args = argparse.ArgumentParser()
    lb2b_args.add_argument("-b", "--bamlist", dest="bamlist", required=True, help="Path to BAM list")
    lb2b_args.add_argument("-r", "--rd3data", dest="rd3data", required=True, help="Path to RD3 data file")
    lb2b_args.add_argument("-k", "--kitlist", dest="kitlist", required=True, help="Path to kit list (kit to bed)")
    lb2b_args.add_argument("-o", "--outfile", dest="outfile", required=True, help="Path to write output file to")
    return vars(lb2b_args.parse_args())


def read_bamlist(bamlistloc):
    bamfiles = {}
    try:
        with open(bamlistloc, 'r') as bamlistfile:
            for fileline in bamlistfile:
                if fileline.strip().endswith((".bam", ".bam.bai", ".bai")):
                    bamsample = fileline.strip().split("/")[-1].split(".")[0]
                    if bamsample not in bamfiles:
                        bamfiles[bamsample] = []
                    bamfiles[bamsample].append(fileline.strip())
    except IOError:
        print("Could not read bam list file")
    finally:
        return bamfiles


def read_rd3data(rd3loc):
    rd3data = {}
    try:
        with open(rd3loc, 'r') as rd3file:
            for fileline in rd3file:
                filelinedata = fileline.strip().split()
                rd3data[filelinedata[1]] = filelinedata[3]
    except IOError:
        print("Could not read rd3 data file")
    finally:
        return rd3data


def read_kitlist(kitlistloc):
    kitfiles = {}
    try:
        with open(kitlistloc, 'r') as kitlistfile:
            for fileline in kitlistfile:
                filelinedata = fileline.strip().split()
                kitfiles[filelinedata[0]] = f"/groups/umcg-solve-rd/tmp01/resources/BED_KITS/{filelinedata[1]}"
    except IOError:
        print("Could not read kit list file")
    finally:
        return kitfiles


def link_bed_to_bam(bamfiles, rd3data, kitfiles, outfileloc):
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("BED\tBAM\n")
            for samplename in rd3data:
                if rd3data[samplename] in kitfiles and samplename in bamfiles:
                    bedfile = kitfiles[rd3data[samplename]]
                    for bamfile in bamfiles[samplename]:
                        outfile.write(f"{bedfile}\t{bamfile}\n")
        file_written = True
    except IOError:
        print("Could not write bed to bam file")
    finally:
        return file_written


def main():
    lb2b_params = get_params()
    print("...READ BAM LIST FILE...")
    bamfiles = read_bamlist(lb2b_params["bamlist"])
    print(f"...READ {len(bamfiles)} BAM LIST ENTRIES")
    print("...READ THE RD3 DATA FILE...")
    rd3data = read_rd3data(lb2b_params["rd3data"])
    print(f"...READ {len(rd3data)} RD3 data entries...")
    print("...READ THE KIT TO BED LIST FILE...")
    kitfiles = read_kitlist(lb2b_params["kitlist"])
    print(f"...READ {len(kitfiles)} KIT ENTRIES...")
    print("...WRITE BED TO BAM FILE...")
    wrote_file = link_bed_to_bam(bamfiles, rd3data, kitfiles, lb2b_params["outfile"])
    print(f"Wrote output file?: {wrote_file}")

if __name__ == "__main__":
    main()
