import argparse
import os


def get_params():
    jobdata_args = argparse.ArgumentParser()
    jobdata_args.add_argument("-d", "--indir", dest="indir", nargs="+", required=True, help="Path to input directory with job output files")
    jobdata_args.add_argument("-o", "--outfile", dest="outfile", required=True, help="Path to write job data output file to")
    return vars(jobdata_args.parse_args())


def gather_job_data(outfiles):
    jobdata = {}
    for outfile in outfiles:
        try:
            with open(outfile, 'r') as joboutfile:
                for fileline in joboutfile:
                    if fileline[0].isdigit():
                        filelinedata = fileline.strip().split()
                        jobdata[filelinedata[0]] = filelinedata
        except IOError:
            print(f"Could not read {outfile}")
    return jobdata


def write_output_file(jobdatastats, outfileloc):
    wrote_file = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("JobID\tElapsed\tAllocCPUS\tAveCPU\tReqMem\tMaxVMSize\tMaxRSS\tMaxDiskRead\tMaxDiskWrite\n")
            for jobstats in jobdatastats:
                outfile.write("\t".join(jobdatastats[jobstats]) +"\n")
        wrote_file = True
    except IOError:
        print(f"Could not write ")
    finally:
        return wrote_file


def get_outfilepaths(indirlocs):
    joboutfiles = []
    for indirloc in indirlocs:
        if not indirloc.endswith("/"):
            indirloc = indirloc + "/"
        indirfiles = os.listdir(indirloc)
        outfiles = [f"{indirloc}{outfile}" for outfile in indirfiles if outfile.endswith(".out")]
        joboutfiles.extend(outfiles)
    return joboutfiles


def main():
    jobdata_params = get_params()
    joboutfiles = get_outfilepaths(jobdata_params["indir"])
    jobstats = gather_job_data(joboutfiles)
    file_written = write_output_file(jobstats, jobdata_params["outfile"])
    print(f"Wrote the output file: {file_written}")


if __name__ == "__main__":
    main()
