import argparse


def get_params():
    denoise_args = argparse.ArgumentParser()
    denoise_args.add_argument("-e", "--eigen-samples", type=int, required=True, dest="eigen-samples", help="Number of eigen samples that should be used")
    denoise_args.add_argument("-j", "--job-out", type=str, required=True, dest="job-outdir", help="Path the job should write output to")
    denoise_args.add_argument("-n", "--name", type=str, required=True, dest="name", help="Name of the batch or bedfile")
    denoise_args.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", help="Path to write script file to")
    denoise_args.add_argument("-s", "--samples-to-pon", type=str, required=True, dest="samples-to-pon", help="Path to file with list of samples to PoN should be used on")
    return vars(denoise_args.parse_args())


def read_sample_list(samplelistfileloc):
    sample_data = {}
    try:
        with open(samplelistfileloc, 'r') as samplelistfile:
            for fileline in samplelistfile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[1] not in sample_data:
                    sample_data[filelinedata[1]] = []
                sample_data[filelinedata[1]].append(filelinedata[0])
    except IOError:
        print("Could not read sample list")
    finally:
        return sample_data


def write_sbatch_script(scriptoutloc, joboutdir, jobname, samplename, samplefile, panelofnormals, eigensamples):
    file_written = False
    try:
        with open(scriptoutloc, 'w') as scriptfile:
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
            scriptfile.write(f"gatk DenoiseReadCounts \\\n")
            scriptfile.write(f"\t-I {samplefile} \\\n")
            scriptfile.write(f"\t--count-panel-of-normals {panelofnormals} \\\n")
            scriptfile.write(f"\t--standardized-copy-ratios {joboutdir}standardized/{samplename}.standardized.tsv \\\n")
            scriptfile.write(f"\t--denoised-copy-ratios {joboutdir}denoised/{samplename}.denoised.tsv \\\n")
            scriptfile.write(f"\t--number-of-eigensamples {eigensamples}\n")
        file_written = True
    except IOError:
        print("Could not write sbatch script file")
    finally:
        return file_written


def write_submitter(submitterfileloc, scriptlist):
    file_written = False
    try:
        with open(submitterfileloc, 'w') as submitterfile:
            for scriptfile in scriptlist:
                submitterfile.write(f"sbatch {scriptfile}\n")
        file_written = True
    except IOError:
        print("Could not write submitter file")
    finally:
        return file_written


def main():
    denoise_params = get_params()

    submitname = denoise_params["name"]
    outdir = denoise_params["outdir"]
    outdir = f"{outdir}/" if not outdir.endswith("/") else outdir
    joboutdir = denoise_params["job-outdir"]
    joboutdir = f"{joboutdir}/" if not joboutdir.endswith("/") else joboutdir
    joboutdir = f"{joboutdir}{submitname}/"
    pon_samples = read_sample_list(denoise_params["samples-to-pon"])

    sbatchlist = []
    for ponfile in pon_samples:
        for ponsample in pon_samples[ponfile]:
            samplename = ponsample.split("/")[-1].split(".")[0]
            scriptoutpath = f"{outdir}drc_{samplename}.sh"
            jobname = f"drc_{samplename}"
            wrote_script = write_sbatch_script(scriptoutpath, joboutdir, jobname, samplename, ponsample, ponfile, denoise_params["eigen-samples"])
            print(f"Wrote job script for {samplename}?: {wrote_script}")
            if wrote_script:
                sbatchlist.append(scriptoutpath)

    submitterpath = f"{outdir}submit_{submitname}.sh"
    wrote_submit = write_submitter(submitterpath, sbatchlist)
    print(f"Wrote submitter script?: {wrote_submit}")


if __name__ == "__main__":
    main()
