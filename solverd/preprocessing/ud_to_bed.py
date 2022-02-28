#!/usr/bin/env python
import os
import argparse


def get_params():
    udtb_args = argparse.ArgumentParser()
    udtb_args.add_argument("-b", "--bed-to-bam", type=str, required=True, dest="bed-to-bam", help="Path to master BED to BAM list file.")
    udtb_args.add_argument("-m", "--master-dir", type=str, dest="master-dir", help="Path to master dir containing BAMs")
    udtb_args.add_argument("-o", "--outdir", type=str, dest="outdir", help="Path to write output file to.")
    udtb_args.add_argument("-r", "--rd-connect", type=str, dest="rd-connect", help="Path to RDConnect file.")
    udtb_args.add_argument("-u", "--ud-dir", type=str, dest="ud-dir", help="Path to directory containing UD samples")
    return vars(udtb_args.parse_args())


def read_bed_to_bam(btbfileloc):
    btbdata = {}
    try:
        with open(btbfileloc, 'r') as btbfile:
            next(btbfile)   # Skip header
            for btbline in btbfile:
                btblinedata = btbline.strip().split("\t")
                if btblinedata[1].endswith(".bam"):
                    bamfile = btblinedata[1].split("/")[-1]
                    bedfile = btblinedata[0].split("/")[-1]
                    btbdata[bamfile] = bedfile
    except IOError:
        print(f"Could not read {btbfileloc}")
    finally:
        return btbdata


def read_rdconnect(rdcfileloc):
    rdcdata = {}
    try:
        with open(rdcfileloc, 'r') as rdcfile:
            for rdcline in rdcfile:
                rdclinedata = rdcline.strip().split("\t")
                rdcdata[rdclinedata[1]] = rdclinedata[0]
    except IOError:
        print("Could not read RDConnect")
    finally:
        return rdcdata


def link_ud_to_bed(bedtobam, udfiles):
    ud_to_bed = {}
    for udfile in udfiles:
        if udfile in bedtobam:
            bedfile = bedtobam[udfile]
            if bedfile not in ud_to_bed:
                ud_to_bed[bedfile] = []
            ud_to_bed[bedfile].append(udfile)
    return ud_to_bed


def write_ud_per_bed(outfileloc, udbedfiles, masterdirloc, uddirloc, kitname):
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            for udfile in udbedfiles:
                outfile.write(f"ln -s {masterdirloc}/{udfile} {uddirloc}/{kitname}/{udfile}\n")
                outfile.write(f"ln -s {masterdirloc}/{udfile}.bai {uddirloc}/{kitname}/{udfile}.bai\n")
        file_written = True
    except IOError:
        print("Could not write ")
    finally:
        return file_written


def main():
    udtb_params = get_params()
    btb_data = read_bed_to_bam(udtb_params["bed-to-bam"])
    print(len(btb_data))
    rdc_data = read_rdconnect(udtb_params["rd-connect"])
    print(len(rdc_data))

    uddir = udtb_params["ud-dir"]
    uddirfiles = os.listdir(uddir)
    ud_samples = [x for x in uddirfiles if x.endswith(".bam")]
    print(len(ud_samples))
    outdir = udtb_params["outdir"]
    masterdir = udtb_params["master-dir"]

    ud_per_bed = link_ud_to_bed(btb_data, ud_samples)
    for bedfile in ud_per_bed:
        kitname = rdc_data[bedfile]
        outfilepath = f"{outdir}/{kitname}.sh"
        write_ud_per_bed(outfilepath, ud_per_bed[bedfile], masterdir, outdir, kitname)


if __name__ == "__main__":
    main()
