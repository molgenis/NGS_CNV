import os
import argparse


def get_params():
    """Define CLI parameters and return the set values.

    Returns
    -------
    dict
        Set parameter values
    """
    udrc_args = argparse.ArgumentParser()
    udrc_args.add_argument("-a", "--awk-script", type=str, required=True, dest="awk-script", help="Path to awk script")
    udrc_args.add_argument("-l", "--intervallist", type=str, required=True, dest="intervallist", help="Path to intervallist to use")
    udrc_args.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", help="Path to directory the generated scripts should be written to")
    udrc_args.add_argument("-s", "--script-outdir", type=str, required=True, dest="script-outdir", help="Path the generated scripts should write their output to")
    udrc_args.add_argument("-u", "--uddir", type=str, required=True, dest="uddir", help="Path to directory with UD BAM files")
    return vars(udrc_args.parse_args())


def read_awk_template(awktemplateloc):
    """Read and return the AWK script that will be used as template code.

    Parameters
    ----------
    awktemplateloc : str
        Path to the AWK script to use as template

    Returns
    -------
    awktemplate : str
        The AWK script
    """
    try:
        awktemplate = ""
        with open(awktemplateloc, 'r') as awkfile:
            awktemplate = awkfile.read()
        return awktemplate
    except IOError:
        print(f"Could not read awk template {awktemplateloc}")


def generate_picard(outfileloc, uddir, udbamfile, intervallistfile, outdir, picardoutfile):
    """Generate Picard CalculateHsMetrics for 

    Parameters
    ----------
    outfileloc : str
        Path to writ this picard script file to
    uddir : str
        Path to directory containing one or more UD BAM files
    udbamfile :
        Path to UD BAM file Picard should use
    intervallistfile : str
        Path to the intervallist Picard should use
    outdir :
        Output directory
    picardoutfile
        Path Picard should write its output file to
    """
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(f"#!/bin/bash\n")
            outfile.write(f"#SBATCH --job-name={udbamfile}_hs\n")
            outfile.write(f"#SBATCH --output={udbamfile}_hs.out\n")
            outfile.write(f"#SBATCH --error={udbamfile}_hs.err\n")
            outfile.write(f"#SBATCH --time=01:00:00\n")
            outfile.write(f"#SBATCH --cpus-per-task=1\n")
            outfile.write(f"#SBATCH --mem=3gb\n")
            outfile.write(f"#SBATCH --nodes=1\n")
            outfile.write(f"#SBATCH --open-mode=append\n")
            outfile.write(f"#SBATCH --export=NONE\n")
            outfile.write(f"#SBATCH --get-user-env=L\n\n")
            outfile.write("module load picard/2.20.5-Java-11-LTS\n")
            outfile.write("module list\n\n")
            outfile.write(f"java -jar -XX:ParallelGCThreads=2 -Xmx2g /apps/software/picard/2.20.5-Java-11-LTS/picard.jar CalculateHsMetrics \\\n")
            outfile.write(f"\tINPUT={uddir}/{udbamfile} \\\n")
            outfile.write(f"\tTARGET_INTERVALS={intervallistfile} \\\n")
            outfile.write(f"\tBAIT_INTERVALS={intervallistfile} \\\n")
            outfile.write(f"\tTMP_DIR={outdir} \\\n")
            outfile.write(f"\tOUTPUT={picardoutfile}\n")
        file_written = True
    except IOError:
        print(f"Could not write Picard script file for {udbamfile}")
    finally:
        return file_written


def generate_awk(outfileloc, awktemplate, udbamfile, awkcsmcfile, awksexfile):
    awkscript = awktemplate
    awkscript = awkscript.replace("${checkSexMeanCoverage}", awkcsmcfile)
    awkscript = awkscript.replace("${dedupBamMetrics}", udbamfile)
    awkscript = awkscript.replace("${whichSex}", awksexfile)

    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(f"#!/bin/bash\n")
            outfile.write(f"#SBATCH --job-name={udbamfile}_awk\n")
            outfile.write(f"#SBATCH --output={udbamfile}_awk.out\n")
            outfile.write(f"#SBATCH --error={udbamfile}_awk.err\n")
            outfile.write(f"#SBATCH --time=00:30:00\n")
            outfile.write(f"#SBATCH --cpus-per-task=1\n")
            outfile.write(f"#SBATCH --mem=2gb\n")
            outfile.write(f"#SBATCH --nodes=1\n")
            outfile.write(f"#SBATCH --open-mode=append\n")
            outfile.write(f"#SBATCH --export=NONE\n")
            outfile.write(f"#SBATCH --get-user-env=L\n\n")
            outfile.write("module list\n\n")
            outfile.write(awkscript)
        file_written = True
    except IOError:
        print("aap")
    finally:
        return file_written


def generate_determine_sex(outfileloc, uddir, udbamfile, intervallistfile, outdir, picardoutfile, awktemplate, awkcsmcfile, awksexfile):
    """Generate Picard CalculateHsMetrics for 

    Parameters
    ----------
    outfileloc : str
        Path to writ this picard script file to
    uddir : str
        Path to directory containing one or more UD BAM files
    udbamfile :
        Path to UD BAM file Picard should use
    intervallistfile : str
        Path to the intervallist Picard should use
    outdir :
        Output directory
    picardoutfile
        Path Picard should write its output file to
    """
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(f"#!/bin/bash\n")
            outfile.write(f"#SBATCH --job-name={udbamfile}_detsex\n")
            outfile.write(f"#SBATCH --output={udbamfile}_detsex.out\n")
            outfile.write(f"#SBATCH --error={udbamfile}_detsex.err\n")
            outfile.write(f"#SBATCH --time=05:30:00\n")
            outfile.write(f"#SBATCH --cpus-per-task=1\n")
            outfile.write(f"#SBATCH --mem=7gb\n")
            outfile.write(f"#SBATCH --nodes=1\n")
            outfile.write(f"#SBATCH --open-mode=append\n")
            outfile.write(f"#SBATCH --export=NONE\n")
            outfile.write(f"#SBATCH --get-user-env=L\n\n")
            outfile.write("module load picard/2.20.5-Java-11-LTS\n")
            outfile.write("module list\n\n")
            outfile.write(f"java -jar -XX:ParallelGCThreads=2 -Xmx2g /apps/software/picard/2.20.5-Java-11-LTS/picard.jar CollectHsMetrics \\\n")
            outfile.write(f"\tINPUT={uddir}/{udbamfile} \\\n")
            outfile.write(f"\tTARGET_INTERVALS={intervallistfile} \\\n")
            outfile.write(f"\tBAIT_INTERVALS={intervallistfile} \\\n")
            outfile.write(f"\tTMP_DIR={outdir} \\\n")
            outfile.write(f"\tOUTPUT={picardoutfile}\n\n\n")

            # Now make the AWK script
            awkscript = awktemplate
            awkscript = awkscript.replace("\"${tmpHsMetricsNonAutosomalRegionChrX}\"", picardoutfile)
            awkscript = awkscript.replace("\"${checkSexMeanCoverage}\"", awkcsmcfile)
            awkscript = awkscript.replace("\"${dedupBamMetrics}\"", udbamfile)
            awkscript = awkscript.replace("\"${whichSex}\"", awksexfile)
            outfile.write(awkscript)
        file_written = True
    except IOError:
        print(f"Could not write script file for {udbamfile}")
    finally:
        return file_written


def main():
    udrc_params = get_params()

    outdir = udrc_params["outdir"]
    intervallist = udrc_params["intervallist"]
    scriptoutdir = udrc_params["script-outdir"]

    uddir = udrc_params["uddir"]
    uddirfiles = os.listdir(uddir)
    udbamfiles = [x for x in uddirfiles if x.endswith(".bam")]

    awk_template = read_awk_template(udrc_params["awk-script"])

    for udbam in udbamfiles:
        picardoutfile = f"{scriptoutdir}/{udbam}.hs_metrics.txt"
        # wrote_picard = generate_picard(f"{outdir}/{udbam}.picard.sh", uddir, udbam, intervallist, outdir, picardoutfile)
        # print(f"Wrote Picard script file for {udbam}?: {wrote_picard}")

        awk_csmc_outfile = f"{scriptoutdir}/{udbam}.csmc.txt"
        awk_sex_outfile = f"{scriptoutdir}/{udbam}.sex.txt"
        # wrote_awk = generate_awk(f"{outdir}/{udbam}.awk.sh", awk_template, )
        # print(f"Wrote awk script file for {udbam}?: {wrote_awk}")

        wrote_script = generate_determine_sex(f"{outdir}/{udbam}.detsex.sh", uddir, udbam, intervallist, outdir, picardoutfile, awk_template, awk_csmc_outfile, awk_sex_outfile)
        print(f"Wrote Picard and AWK script file for {udbam}?: {wrote_script}")


if __name__ == "__main__":
    main()
