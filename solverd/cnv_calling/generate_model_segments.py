#!/usr/bin/env python
import os
import argparse


def get_params():
    ms_args = argparse.ArgumentParser()
    ms_args.add_argument("-c", "--cac-dir", type=str, dest="cac-dir", required=True, help="Path to directory")
    ms_args.add_argument("-d", "--drc-dir", type=str, dest="drc-dir", required=True, help="")
    ms_args.add_argument("-o", "--outdir", type=str, dest="outdir", required=True, help="Path to write script files to")
    ms_args.add_argument("-j", "--job-outdir", type=str, dest="job-outdir", required=True, help="Path the job should write the output to")
    # ms_args.add_argument("-n", "--name", type=str, dest="name", required=True, help="Batch name")
    return vars(ms_args.parse_args())


def write_job_script(scriptoutpath, joboutpath, joboutprefix, jobname, drcfile, cacfile):
    file_written = False
    try:
        with open(scriptoutpath, 'w') as scriptfile:
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

            # Write the ModelSegments command and the required parameters
            scriptfile.write("gatk ModelSegments \\\n")
            scriptfile.write(f"\t--denoised-copy-ratios {drcfile} \\\n")
            scriptfile.write(f"\t--allelic-counts {cacfile} \\\n")
            scriptfile.write(f"\t--output {joboutpath} \\\n")
            scriptfile.write(f"\t--output-prefix {joboutprefix}\n")
        file_written = True
    except IOError:
        print("Could not write script file")
    finally:
        return file_written


def write_submitter(submitscriptloc, jobscripts):
    file_written = False
    try:
        with open(submitscriptloc, 'w') as submitfile:
            for jscript in jobscripts:
                submitfile.write(f"sbatch {jscript}\n")
        file_written = True
    except IOError:
        print("Could not write job submitter script")
    finally:
        return file_written


def get_input_files(inputfilesdir, fileext):
    indirfiles = os.listdir(inputfilesdir)
    infiles = {x.split(".")[0]:f"{inputfilesdir}{x}" for x in indirfiles if x.endswith(fileext)}
    return infiles


def main():
    ms_params = get_params()
    drcdir = ms_params["drc-dir"]
    drcdir = f"{drcdir}/" if not drcdir.endswith("/") else drcdir
    cacdir = ms_params["cac-dir"]
    cacdir = f"{cacdir}/" if not cacdir.endswith("/") else cacdir
    outdir = ms_params["outdir"]
    outdir = f"{outdir}/" if not outdir.endswith("/") else outdir
    joboutdir = ms_params["job-outdir"]
    joboutdir = f"{joboutdir}/" if not joboutdir.endswith("/") else outdir
    # batchname = ms_params["name"]

    drc_files = get_input_files(drcdir, ".tsv")
    print(f"Found {len(drc_files)} denoised files")
    cac_files = get_input_files(cacdir, ".tsv")
    print(f"Found {len(cac_files)} allelic count files")
    valid_samples = set(drc_files.keys()) & set(cac_files.keys())
    invalid_samples = set(drc_files.keys()) ^ set(cac_files.keys())

    # Start making the job scripts
    submitter_list = []
    for samplename in valid_samples:
        scriptfilepath = f"{outdir}ms_{samplename}.sh"
        jobname = f"ms_{samplename}"
        wrote_file = write_job_script(scriptfilepath, joboutdir, samplename, jobname, drc_files[samplename], cac_files[samplename])
        print(f"Wrote job script {scriptfilepath}?: {wrote_file}")

        if wrote_file:
            submitter_list.append(scriptfilepath)

    submitfileloc = f"{outdir}submit_jobs.sh"
    wrote_submitter = write_submitter(submitfileloc, submitter_list)
    print(f"Wrote submitter file?: {wrote_submitter}")

    # Check if there are any invalid samples for which no job script could be written
    if len(invalid_samples) > 0:
        print(f"Could not write job scripts for samples: {invalid_samples}")


if __name__ == "__main__":
    main()
