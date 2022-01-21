import os
import argparse


def get_params():
    pms_args = argparse.ArgumentParser()
    pms_args.add_argument("-a", "--allelic-dir", type=str, dest="allelic-dir", required=True, help="Path to directory containing the allelic counts files")
    pms_args.add_argument("-d", "--denoised-dir", type=str, dest="denoised-dir", required=True, help="Path to directory containing the denoised read counts")
    pms_args.add_argument("-j", "--joboutdir", type=str, dest="joboutdir", required=True, help="Path where the job should write the output to")
    pms_args.add_argument("-m", "--modelsegments-dir", type=str, dest="modelsegments-dir", required=True, help="Path to directory containing model segments files")
    pms_args.add_argument("-o", "--outdir", type=str, dest="outdir", required=True, help="Path to write script to")
    pms_args.add_argument("-s", "--sequence-dict", type=str, dest="sequence-dict", required=True, help="Path to sequence directory to use")
    return vars(pms_args.parse_args())


def write_job_script(scriptfileloc, drcfile, cacfile, msfile, seqdict, joboutdir, jobname):
    file_written = False
    try:
        with open(scriptfileloc, 'w') as scriptfile:
            # Write the header
            scriptfile.write("#!/bin/bash\n")
            scriptfile.write(f"#SBATCH --job-name={jobname}\n")
            scriptfile.write(f"#SBATCH --output={jobname}.out\n")
            scriptfile.write(f"#SBATCH --error={jobname}.err\n")
            scriptfile.write("#SBATCH --time=00:20:00\n")
            scriptfile.write("#SBATCH --cpus-per-task=1\n")
            scriptfile.write("#SBATCH --mem=5gb\n")
            scriptfile.write("#SBATCH --nodes=1\n")
            scriptfile.write("#SBATCH --open-mode=append\n")
            scriptfile.write("#SBATCH --export=NONE\n")
            scriptfile.write("#SBATCH --get-user-env=L\n\n")
            scriptfile.write("module load GATK/4.1.4.1-Java-8-LTS\n")
            scriptfile.write("module load RPlus/4.0.3-foss-2018b-v21.08.1\n")
            scriptfile.write("module list\n\n")

            # Write the command
            scriptfile.write("gatk PlotModeledSegments \\\n")
            scriptfile.write(f"--denoised-copy-ratios {drcfile} \\\n")
            # scriptfile.write(f"--allelic-counts {cacfile} \\\n")
            scriptfile.write(f"--segments {msfile} \\\n")
            scriptfile.write(f"--sequence-dictionary {seqdict} \\\n")
            scriptfile.write(f"--output-prefix {jobname} \\\n")
            scriptfile.write(f"--output {joboutdir} \n")
        file_written = True
    except IOError:
        print(f"Could not write job script for {scriptfileloc}")
    finally:
        return file_written


def write_submitter(submitscriptloc, jobscripts):
    file_written = False
    try:
        with open(submitscriptloc, 'w') as submitscript:
            for jobscript in jobscripts:
                submitscript.write(f"sbatch {jobscript}\n")
        file_written = True
    except IOError:
        print("Could not write submitter script")
    finally:
        return file_written


def main():
    pms_params = get_params()

    # Gather the directory paths
    cacdir = pms_params["allelic-dir"]
    drcdir = pms_params["denoised-dir"]
    joboutdir = pms_params["joboutdir"]
    msdir = pms_params["modelsegments-dir"]
    outdir = pms_params["outdir"]

    # Add a trailing / to each directory path
    cacdir = f"{cacdir}/" if not cacdir.endswith("/") else cacdir
    drcdir = f"{drcdir}/" if not drcdir.endswith("/") else drcdir
    joboutdir = f"{joboutdir}/" if not joboutdir.endswith("/") else joboutdir
    msdir = f"{msdir}/" if not msdir.endswith("/") else msdir
    outdir = f"{outdir}/" if not outdir.endswith("/") else outdir

    # Collect the input files from the cac, drc and ms directories
    cac_files = {x.split(".")[0]:x for x in os.listdir(cacdir) if x.endswith(".cac.tsv")}
    drc_files = {x.split(".")[0]:x for x in os.listdir(drcdir) if x.endswith(".denoised.tsv")}
    ms_files = {x.split(".")[0]:x for x in os.listdir(msdir) if x.endswith(".modelFinal.seg")}

    # Determine the valid and invalid samples
    valid_samples = set(cac_files.keys()) & set(drc_files.keys()) & set(ms_files.keys())
    invalid_cac_samples = set(cac_files.keys()) - valid_samples
    invalid_drc_samples = set(drc_files.keys()) - valid_samples
    invalid_ms_samples = set(ms_files.keys()) - valid_samples

    # Start making the job scripts
    written_jobscripts = []
    for samplename in valid_samples:
        jobname = f"pms_{samplename}"
        scriptoutloc = f"{outdir}{jobname}.sh"
        wrote_file = write_job_script(scriptoutloc, f"{drcdir}{drc_files[samplename]}", f"{cacdir}{cac_files[samplename]}", f"{msdir}{ms_files[samplename]}", pms_params["sequence-dict"], joboutdir, jobname)
        print(f"Wrote jobscript {scriptoutloc}?: {wrote_file}")

        if wrote_file:
            written_jobscripts.append(scriptoutloc)

    # Make the submitter script
    if len(written_jobscripts) > 0:
        wrote_submitter = write_submitter(f"{outdir}submit_jobs.sh", written_jobscripts)
        print(f"Wrote submitter script?: {wrote_submitter}")

    # Print out the invalid samples if there are any
    if len(invalid_cac_samples) > 0:
        print(f"Invalid CAC samples: {invalid_cac_samples}")

    if len(invalid_drc_samples) > 0:
        print(f"Invalid DRC samples: {invalid_drc_samples}")

    if len(invalid_ms_samples) > 0:
        print(f"Invalid MS samples: {invalid_ms_samples}")


if __name__ == "__main__":
    main()
