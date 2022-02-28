#!/usr/bin/env python
import os
import argparse


def get_params():
    pon_args = argparse.ArgumentParser()
    pon_args.add_argument("-b", "--bam-to-sex", type=str, required=True, dest="bam-to-sex", help="Path to BAM to sex file")
    pon_args.add_argument("-c", "--s4-clusters", type=str, required=True, dest="s4-clusters", help="Path to combined S4 clusters file")
    pon_args.add_argument("-j", "--job-outdir", type=str, required=True, dest="job-outdir", help="Where job should write output to")
    pon_args.add_argument("-n", "--name", type=str, required=True, dest="name", help="Name to use for output")
    pon_args.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", help="Path to output directory")
    pon_args.add_argument("-r", "--read-counts-dir", type=str, required=True, dest="read-counts-dir", help="Path to directory with read counts files")
    pon_args.add_argument("-s", "--script-oudir", type=str, required=True, dest="script-outdir", help="Path to write scripts to")
    return vars(pon_args.parse_args())


def read_bam_to_sex(btsfileloc):
    btsdata = {"F": [], "M": []}
    try:
        with open(btsfileloc, 'r') as btsfile:
            next(btsfile)
            for btsline in btsfile:
                btslinedata = btsline.strip().split("\t")
                if len(btslinedata) == 3:
                    if btslinedata[2] == "F":
                        btsdata["F"].append(btslinedata[0])
                    if btslinedata[2] == "M":
                        btsdata["M"].append(btslinedata[0])
    except IOError:
        print(f"Could not open")
    finally:
        return btsdata


def read_combineds4_data(combineds4_fileloc):
    s4data = {}
    try:
        with open(combineds4_fileloc, 'r') as cs4file:
            for cs4line in cs4file:
                if cs4line.startswith("["):
                    s4data = get_cluster_samples(cs4line.strip()[1:-1], s4data)
    except IOError:
        print("Could not read S4 combined cluster data")
    finally:
        return s4data


def get_cluster_samples(clusterdata, batchclusters):
    cluster_samples = clusterdata.split("),")
    for csample in cluster_samples:
        clusternum = csample.strip().split(",")[1].strip()

        if clusternum.endswith(")"):
            clusternum = clusternum[0:-1]

        if clusternum not in batchclusters:
            batchclusters[clusternum] = []
        samplename = csample.strip().split(",")[0][2:-1]
        batchclusters[clusternum].append(samplename)
    return batchclusters


def fetch_fm_samples(ccluster, bamtosex):
    """Fetch all male/female samples for a S4 ClusterWES cluster.

    Parameters
    ----------
    ccluster : list of str
        List of cluster samples
    bamtosex : list of str
        Liset of all male/female samples
    """
    fmsamples = []
    for clustersample in ccluster:
        if clustersample in bamtosex:
            fmsamples.append(clustersample)
    return fmsamples


def get_readcount_files(readcountsdir):
    """Collect and return the

    Parameters
    ----------
    readcountsdir : str
        Path to directory containing the HDF5 read count files

    Returns
    -------
    """
    readcountsdir = f"{readcountsdir}/" if not readcountsdir.endswith("/") else readcountsdir
    rcdirfiles = os.listdir(readcountsdir)
    hdf5files = [f"{readcountsdir}{hdf5file}" for hdf5file in rcdirfiles]
    return hdf5files


def write_pon_job(ponsamples, readcountfiles, ponscriptoutloc, ponfileoutloc, jobname):
    """Write a job script for creating the PanelOfNormals

    Parameters
    ----------
    ponsamples : list of str
        Samples to use for making the PanelOfNormals
    readcountfiles : list of str
        Read count files
    ponscriptoutloc : str
        Path to write the PoN job script to
    ponfileoutloc : str
        Path job should write output file to
    jobname : str
        Name for the PoN job

    Returns
    -------
    file_written : bool
        True if file has been successfully written, False if not
    """
    file_written = False
    try:
        with open(ponscriptoutloc, 'w') as ponjobfile:
            ponjobfile.write("#!/bin/bash\n")
            ponjobfile.write(f"#SBATCH --job-name={jobname}\n")
            ponjobfile.write(f"#SBATCH --output={jobname}.out\n")
            ponjobfile.write(f"#SBATCH --error={jobname}.err\n")
            ponjobfile.write("#SBATCH --time=01:00:00\n")
            ponjobfile.write("#SBATCH --cpus-per-task=1\n")
            ponjobfile.write("#SBATCH --mem=5gb\n")
            ponjobfile.write("#SBATCH --nodes=1\n")
            ponjobfile.write("#SBATCH --open-mode=append\n")
            ponjobfile.write("#SBATCH --export=NONE\n")
            ponjobfile.write("#SBATCH --get-user-env=L\n\n")
            ponjobfile.write("module load GATK/4.1.4.1-Java-8-LTS\n")
            ponjobfile.write("module list\n\n")
            ponjobfile.write("gatk CreateReadCountPanelOfNormals \\\n")

            for rcfile in readcountfiles:
                samplename = rcfile.split("/")[-1].split(".")[0]
                if samplename in ponsamples:
                    ponjobfile.write(f"\t-I {rcfile} \\\n")
            ponjobfile.write(f"\t-O {ponfileoutloc}\n")
        file_written = True
    except IOError:
        print(f"")
    finally:
        return file_written


def write_samples_to_pon(ponsamples, readcountfiles, ponfile, outfileloc):
    """Write which PanelOfNormals to use for a set of samples.

    Parameters
    ----------
    ponsamples : list of str
        Samples used to construct the PanelOfNormals
    ponfile : str
        Path to the PanelOfNormals file
    outfileloc : str
        Path to write output file to

    Returns
    -------
    file_written : bool
        True if output file has successfully been written, False if not
    """
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            for rcfile in readcountfiles:
                samplename = rcfile.split("/")[-1].split(".")[0]
                if samplename in ponsamples:
                    outfile.write(f"{rcfile}\t{ponfile}\n")
        file_written = True
    except IOError:
        print(f"Could not write output file {outfileloc}")
    finally:
        return file_written


def main():
    pon_params = get_params()

    # Make sure that paths to directories end in a /
    prefix_name = pon_params["name"]
    outdir = pon_params["outdir"]
    outdir = f"{outdir}/" if not outdir.endswith("/") else outdir
    joboutdir = pon_params["job-outdir"]
    joboutdir = f"{joboutdir}/" if not joboutdir.endswith("/") else joboutdir
    readcountsdir = pon_params["read-counts-dir"]
    readcountsdir = f"{readcountsdir}/" if not readcountsdir.endswith("/") else readcountsdir
    scriptoutdir = pon_params["script-outdir"]
    scriptoutdir = f"{scriptoutdir}/" if not scriptoutdir.endswith("/") else scriptoutdir

    # Read the input data files
    btsdata = read_bam_to_sex(pon_params["bam-to-sex"])
    s4cdata = read_combineds4_data(pon_params["s4-clusters"])
    readcountfiles = get_readcount_files(readcountsdir)

    # Iterate over the combined clusters
    for clusternum in s4cdata:
        # Make some outfile path
        pon_name_f = f"{joboutdir}{prefix_name}_fpon_{clusternum}.hdf5"
        pon_name_m = f"{joboutdir}{prefix_name}_mpon_{clusternum}.hdf5"
        ponjob_name_f = f"{prefix_name}_fpon_{clusternum}"
        ponjob_name_m = f"{prefix_name}_mpon_{clusternum}"
        script_path_f = f"{scriptoutdir}{prefix_name}_makefpon_{clusternum}.sh"
        script_path_m = f"{scriptoutdir}{prefix_name}_makempon_{clusternum}.sh"
        ponlink_fpath = f"{outdir}{prefix_name}_{clusternum}_fsamples.txt"
        ponlink_mpath = f"{outdir}{prefix_name}_{clusternum}_msamples.txt"

        # Select the male and female samples for the PanelOfNormals
        fsamples = fetch_fm_samples(s4cdata[clusternum], btsdata["F"])
        msamples = fetch_fm_samples(s4cdata[clusternum], btsdata["M"])

        # Make the PanelOfNormals job with the samples
        wrote_ponjob_f = write_pon_job(fsamples, readcountfiles, script_path_f, pon_name_f, ponjob_name_f)
        print(f"Successfully wrote female PoN jobscript for cluster {clusternum}?: {wrote_ponjob_f}")
        wrote_ponjob_m = write_pon_job(msamples, readcountfiles, script_path_m, pon_name_m, ponjob_name_m)
        print(f"Successfully wrote male PoN jobscript for cluster {clusternum}?: {wrote_ponjob_m}")

        # Make the file listing the PanelOfNormals to use for each read counts file
        wrote_ponlink_f = write_samples_to_pon(fsamples, readcountfiles, pon_name_f, ponlink_fpath)
        print(f"Successfully wrote female sample to PoN file for cluster {clusternum}?: {wrote_ponlink_f}")
        wrote_ponlink_m = write_samples_to_pon(msamples, readcountfiles, pon_name_m, ponlink_mpath)
        print(f"Successfully wrote male sample to PoN file for cluster {clusternum}?: {wrote_ponlink_m}")


if __name__ == "__main__":
    main()

# USAGE:
# python /groups/umcg-solve-rd/tmp01/umcg-mbeukers/solverd_scripts/solve_rd_make_pon.py
# -b /groups/umcg-solve-rd/tmp01/umcg-mbeukers/rd3_data/sex_to_bam.txt
# -c /groups/umcg-solve-rd/tmp01/umcg-mbeukers/combined_s4_clusterwes/mergesplit/batch3.txt
# -j /groups/umcg-solve-rd/tmp01/umcg-mbeukers/gatk4_freeze1/pon/batch3
# -n batch3
# -o /groups/umcg-solve-rd/tmp01/umcg-mbeukers/gatk4_freeze1/data2/samples_to_pon
# -r /groups/umcg-solve-rd/tmp01/umcg-mbeukers/gatk4_freeze1/crc/batch3
# -s /groups/umcg-solve-rd/tmp01/umcg-mbeukers/gatk4_freeze1/scripts/jobs/pon_jobs/batch3
