import argparse
import os


def get_params():
    pdcr_args = argparse.ArgumentParser()
    pdcr_args.add_argument("-d", "--drcdir", type=str, dest="drcdir", required=True, help="Path to directory containing denoised and standardized data in two subfolders")
    pdcr_args.add_argument("-j", "--joboutdir", type=str, dest="joboutdir", required=True, help="Path where each job should write its output to")
    pdcr_args.add_argument("-o", "--outdir", type=str, dest="outdir", required=True, help="Path to write script files to")
    pdcr_args.add_argument("-s", "--sequence-dictionary", type=str, dest="sequence-dictionary", required=True, help="")
    return vars(pdcr_args.parse_args())


def write_job_script(scriptfileloc, drcdir, denoisedfile, standardizedfile, seqdict, jobname, joboutdir, joboutprefix):
    file_written = False
    try:
        with open(scriptfileloc, 'w') as scriptfile:
            # Write the SBATCH header
            scriptfile.write("#!/bin/bash\n")
            scriptfile.write(f"#SBATCH --job-name={jobname}\n")
            scriptfile.write(f"#SBATCH --output={jobname}.out\n")
            scriptfile.write(f"#SBATCH --error={jobname}.err\n")
            scriptfile.write("#SBATCH --time=00:30:00\n")
            scriptfile.write("#SBATCH --cpus-per-task=1\n")
            scriptfile.write("#SBATCH --mem=5gb\n")
            scriptfile.write("#SBATCH --nodes=1\n")
            scriptfile.write("#SBATCH --open-mode=append\n")
            scriptfile.write("#SBATCH --export=NONE\n")
            scriptfile.write("#SBATCH --get-user-env=L\n\n")
            scriptfile.write("module load GATK/4.1.4.1-Java-8-LTS\n")
            scriptfile.write("module load RPlus/4.0.3-foss-2018b-v21.08.1\n")
            scriptfile.write("module list\n\n")

            # Write the job command
            scriptfile.write("gatk PlotDenoisedCopyRatios \\\n")
            scriptfile.write(f"--denoised-copy-ratios {drcdir}denoised/{denoisedfile} \\\n")
            scriptfile.write(f"--standardized-copy-ratios {drcdir}standardized/{standardizedfile} \\\n")
            scriptfile.write(f"--sequence-dictionary {seqdict} \\\n")
            scriptfile.write(f"--output {joboutdir} \\\n")
            scriptfile.write(f"--output-prefix {joboutprefix}\n")
        file_written = True
    except IOError:
        print("Could not write job script {scriptifleloc}")
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
        print("Could not write submitter file")
    finally:
        return file_written


def main():
    pdcr_params = get_params()

    # Make directory paths
    drcdir = pdcr_params["drcdir"]
    drcdir = f"{drcdir}/" if not drcdir.endswith("/") else drcdir
    joboutdir = pdcr_params["joboutdir"]
    joboutdir = f"{joboutdir}/" if not joboutdir.endswith("/") else joboutdir
    outdir = pdcr_params["outdir"]
    outdir = f"{outdir}/" if not outdir.endswith("/") else outdir

    # Gather denosied and standardized files
    denoised_files = {x.split(".")[0]:x for x in os.listdir(f"{drcdir}denoised/") if x.endswith(".denoised.tsv")}
    standardized_files = {x.split(".")[0]:x for x in os.listdir(f"{drcdir}standardized/") if x.endswith(".standardized.tsv")}

    # Start writing job scripts, one for each sample
    scripts_written = []
    samples_to_process = set(denoised_files.keys()) & set(standardized_files.keys())
    for samplename in samples_to_process:
        jobname = f"pdcr_{samplename}"
        scriptfilepath = f"{outdir}{jobname}.sh"
        wrote_script = write_job_script(scriptfilepath, drcdir, denoised_files[samplename], standardized_files[samplename], pdcr_params["sequence-dictionary"], jobname, joboutdir, samplename)
        print(f"Wrote job script for {samplename}?: {wrote_script}")

        if wrote_script:
            scripts_written.append(scriptfilepath)

    # Write the submitter script.
    if len(scripts_written) > 0 :
        wrote_submitter = write_submitter(f"{outdir}submit_jobs.sh", scripts_written)
        print(f"Wrote job submitter script?: {wrote_submitter}")


if __name__ == "__main__":
    main()
