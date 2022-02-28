#!/usr/bin/env python
import os
import argparse


def get_params():
    ccrs_args = argparse.ArgumentParser()
    ccrs_args.add_argument("-m", "--modelsegments-dir", type=str, dest="modelsegments-dir", required=True, help="")
    ccrs_args.add_argument("-j", "--joboutdir", type=str, dest="joboutdir", required=True, help="")
    ccrs_args.add_argument("-o", "--outdir", type=str, dest="outdir", required=True, help="")
    return vars(ccrs_args.parse_args())


def write_job_script(scriptfilepath, msfile, ccrsfile, jobname):
    file_written = False
    try:
        with open(scriptfilepath, 'w') as scriptfile:
            # Write the SBATCH header
            scriptfile.write("#!/bin/bash\n")
            scriptfile.write(f"#SBATCH --job-name={jobname}\n")
            scriptfile.write(f"#SBATCH --output={jobname}.out\n")
            scriptfile.write(f"#SBATCH --error={jobname}.err\n")
            scriptfile.write("#SBATCH --time=00:15:00\n")
            scriptfile.write("#SBATCH --cpus-per-task=1\n")
            scriptfile.write("#SBATCH --mem=5gb\n")
            scriptfile.write("#SBATCH --nodes=1\n")
            scriptfile.write("#SBATCH --open-mode=append\n")
            scriptfile.write("#SBATCH --export=NONE\n")
            scriptfile.write("#SBATCH --get-user-env=L\n\n")
            scriptfile.write("module load GATK/4.1.4.1-Java-8-LTS\n")
            scriptfile.write("module list\n\n")

            # Write the CallCopyatioSegments command with required parameters
            scriptfile.write("gatk CallCopyRatioSegments \\\n")
            scriptfile.write(f"--input {msfile} \\\n")
            scriptfile.write(f"--output {ccrsfile}\n")
        file_written = True
    except IOError:
        print(f"Could not write jobscript {scriptfilepath}")
    finally:
        return file_written


def write_submitter(submitterpath, jobscripts):
    file_written = False
    try:
        with open(submitterpath, 'w') as submitfile:
            for jobscript in jobscripts:
                submitfile.write(f"sbatch {jobscript}\n")
        file_written = True
    except IOError:
        print("Could not write submit script")
    finally:
        return file_written


def main():
    ccrs_params = get_params()

    msdir = ccrs_params["modelsegments-dir"]
    joboutdir = ccrs_params["joboutdir"]
    outdir = ccrs_params["outdir"]

    msdir = f"{msdir}/" if not msdir.endswith("/") else msdir
    joboutdir = f"{joboutdir}/" if not joboutdir.endswith("/") else joboutdir
    outdir = f"{outdir}/" if not outdir.endswith("/") else outdir

    ms_files = [x for x in os.listdir(msdir) if x.endswith(".cr.seg")]

    submitter_list = []
    for msfile in ms_files:
        samplename = msfile.split(".")[0]
        jobscriptpath = f"{outdir}ccrs_{samplename}.sh"
        wrote_file = write_job_script(jobscriptpath, f"{msdir}{msfile}", f"{joboutdir}{samplename}.called.seg", f"ccrs_{samplename}")
        print(f"Wrote jobscript {jobscriptpath}?: {wrote_file}")

        if wrote_file:
            submitter_list.append(jobscriptpath)

    if len(submitter_list) > 0:
        wrote_file = write_submitter(f"{outdir}submitter_jobs.sh", submitter_list)
        print(f"Wrote submit script?: {wrote_file}")


if __name__ == "__main__":
    main()
