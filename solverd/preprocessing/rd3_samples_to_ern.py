#!/usr/bin/env python
import sys


def read_rd3_exp(expfileloc):
    expdata = {}
    try:
        with open(expfileloc, 'r') as expfile:
            next(expfile)
            for fileline in expfile:
                filelinedata = fileline.strip().split("\t")
                expdata[filelinedata[0]] = filelinedata[1]
    except IOError:
        print("Could not read RD3 experiment data")
    finally:
        return expdata


def read_rd3_sam(samfileloc):
    samdata = {}
    try:
        with open(samfileloc, 'r') as samfile:
            next(samfile)
            for fileline in samfile:
                filelinedata = fileline.strip().split("\t")
                if len(filelinedata) == 2:
                    samdata[filelinedata[0]] = filelinedata[1]
    except IOError:
        print("Could not read RD3 samples data")
    finally:
        return samdata


def read_rd3_sub(subfileloc):
    subdata = {}
    try:
        with open(subfileloc, 'r') as subfile:
            next(subfile)
            for fileline in subfile:
                filelinedata = fileline.strip().split("\t")
                subdata[filelinedata[0]] = filelinedata[-1]
    except IOError:
        print("Could not read RD3 subjects data")
    finally:
        return subdata


def write_sample_to_ern(expdata, samdata, subdata, outfileloc):
    try:
        with open(outfileloc, 'w') as outfile:
            for ename in expdata:
                if expdata[ename] in samdata:
                    if samdata[expdata[ename]] in subdata:
                        outfile.write(f"{ename}\t{subdata[samdata[expdata[ename]]]}\n")
    except IOError:
        print("Could not write sample to ERN file")


# Read the data
exp_data = read_rd3_exp(sys.argv[1])
sam_data = read_rd3_sam(sys.argv[2])
sub_data = read_rd3_sub(sys.argv[3])
write_sample_to_ern(exp_data, sam_data, sub_data, sys.argv[4])
