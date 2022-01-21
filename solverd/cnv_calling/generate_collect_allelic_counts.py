import os
import argparse


def get_params():
    cac_args = argparse.ArgumentParser()
    cac_args.add_argument("-i", "--indir", type=str, dest="indir", required=True, help="Path to the directory with BAM files for a single BED file")
    cac_args.add_argument("-j", "--job-outdir", type=str, required=True, dest="job-outdir", help="Path to folder the job should write the output to")
    cac_args.add_argument("-o", "--outdir", type=str, dest="outdir", help="Path to directory to write script files to")
    cac_args.add_argument("-p", "--prefix", type=str, required=True, dest="prefix", help="Prefix to use for job and output")
    cac_args.add_argument("-r", "--genome-reference", type=str, required=True, dest="genome-reference", help="Path to genome reference to use")
    cac_args.add_argument("-s", "--snp-data", type=str, required=True, dest="snp-data", help="Path to gnomad snp data to use")
    return vars(cac_args.parse_args())


def write_job_script(scriptoutpath, joboutpath, jobname, bamfile, genomereference, snpdata):
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

            # Write the CollectAllelicCounts command and the required parameters
            scriptfile.write("gatk CollectAllelicCounts \\\n")
            scriptfile.write(f"\t-L {snpdata} \\\n")
            scriptfile.write(f"\t-I {bamfile} \\\n")
            scriptfile.write(f"\t-R {genomereference} \\\n")
            scriptfile.write(f"\t-O {joboutpath} \\\n")

            # Write the read filter parameters
            scriptfile.write("\t--read-filter NonZeroReferenceLengthAlignmentReadFilter \\\n")
            scriptfile.write("\t--read-filter MappedReadFilter \\\n")
            scriptfile.write("\t--read-filter NotDuplicateReadFilter \\\n")
            scriptfile.write("\t--read-filter WellformedReadFilter \\\n")
            scriptfile.write("\t--read-filter ProperlyPairedReadFilter \n")
        file_written = True
    except IOError:
        print(f"Could not write script file {scriptoutpath}")
    finally:
        return file_written


def write_submitter(submitterloc, job_scripts):
    file_written = False
    try:
        with open(submitterloc, 'w') as submitfile:
            for jscript in job_scripts:
                submitfile.write(f"sbatch {jscript}\n")
        file_written = True
    except IOError:
        print(f"Could not write submitter script {submitterloc}")
    finally:
        return file_written


def main():
    cac_params = get_params()

    # Ensure a trailing / to directory paths
    bamdir = cac_params["indir"]
    bamdir = f"{bamdir}/" if not bamdir.endswith("/") else bamdir
    outdir = cac_params["outdir"]
    outdir = f"{outdir}/" if not outdir.endswith("/") else outdir
    joboutdir = cac_params["job-outdir"]
    joboutdir = f"{joboutdir}/" if not joboutdir.endswith("/") else joboutdir

    # Gather the BAM files
    bamfiles = [f"{bamdir}{bamfile}" for bamfile in os.listdir(bamdir) if bamfile.endswith(".bam")]

    # Start writing the job scripts
    cac_job_scripts = []
    for bamfile in bamfiles:
        samplename = bamfile.split("/")[-1].split(".")[0]
        outprefix = cac_params["prefix"]
        jobname = f"{outprefix}_{samplename}"
        script_out_path = f"{outdir}{samplename}_cac.sh"
        job_out_path = f"{joboutdir}{samplename}.cac.tsv"

        wrote_job_script = write_job_script(script_out_path, job_out_path, jobname, bamfile, cac_params["genome-reference"], cac_params["snp-data"])
        if wrote_job_script:
            print(f"Wrote job script {script_out_path}?: {wrote_job_script}")
            cac_job_scripts.append(script_out_path)

    # Write the submitter script
    submitfile = f"{outdir}submit_{outprefix}.sh"
    wrote_submitter = write_submitter(submitfile, cac_job_scripts)
    print(f"Wrote job submitter script?: {wrote_submitter}")


if __name__ == "__main__":
    main()
