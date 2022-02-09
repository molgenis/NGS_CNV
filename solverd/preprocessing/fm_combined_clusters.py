import argparse


def get_params():
    """Define CLI parameter and return the set values.

    Returns
    -------
    dict
        Set CLI parameters values
    """
    fmcc_args = argparse.ArgumentParser()
    fmcc_args.add_argument("-b", "--bam-to-sex", type=str, required=True, dest="bam-to-sex", help="Path to BAM to sex file made bij link_sex_to_bam.py")
    fmcc_args.add_argument("-c", "--combined-s4", type=str, required=True, dest="combined-s4", help="Path to combined S4 ClusterWES file")
    fmcc_args.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", help="Path to write output files to")
    fmcc_args.add_argument("-p", "--pon-out", type=str, required=True, dest="pon-out", help="Path where the PoN should be written to by GATK4")
    fmcc_args.add_argument("-r", "--read-counts-dir", type=str, required=True, dest="read-counts-dir", help="Path to directory with read count files")
    return vars(fmcc_args.parse_args())


def read_bamsex_data(bamsex_fileloc):
    """Read and return BAM to sex data.

    Parameters
    ----------
    bamsex_fileloc : str
        Path to BAM to sex file

    Returns
    -------
    bamsexdata : dict
        BAM names and file paths for F and M
    """
    bamsexdata = {"F": {}, "M": {}}
    try:
        with open(bamsex_fileloc, 'r') as bsfile:
            for bsline in bsfile:
                bslinedata = bsline.strip().split("\t")
                if bslinedata[2] == "F":
                    bamsexdata["F"][bslinedata[0]] = bslinedata[1]
                elif bslinedata[2] == "M":
                    bamsexdata["M"][bslinedata[0]] = bslinedata[1]
    except IOError:
        print("aap")
    finally:
        return bamsexdata


def read_bamsex_data_2(bamsex_fileloc):
    bamsexdata = {"F": [], "M": []}
    try:
        with open(bamsex_fileloc, 'r') as bsfile:
            for bsline in bsfile:
                bslinedata = bsline.strip().split("\t")
                if len(bslinedata) == 3:
                    if bslinedata[1] != "NA":
                        if bslinedata[2] == "F":
                            bamsexdata["F"].append(bslinedata[0])
                        elif bslinedata[2] == "M":
                            bamsexdata["M"].append(bslinedata[0])
    except IOError:
        print("Could not read BAM to sex file")
    finally:
        return bamsexdata


def read_combineds4_data(combineds4_fileloc):
    """Read and return combined S4 ClusterWES data.

    Parameters
    ----------
    combineds4_fileloc : str
        Path to combined S4 ClusterWES data

    Returns
    -------
    """
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


def make_fm_pon(clustersamples, bamsexdata):
    """Places male and female samples in different panel of normals for a single cluster

    Parameters
    ----------
    clustersamples : list of str
        Sample names in a cluster
    bamsexdata : dict
        BAM files sorted per sex

    Returns
    -------
    cluster_fmpons : dict
        List of S4 cluster samples for female and male panel of normals
    """
    cluster_fmpons = {"F": [], "M": []}

    # Itereate over the female samples
    for fsample in bamsexdata["F"]:
        if fsample in clustersamples:
            cluster_fmpons["F"].append(fsample)

    # Iterate over the male samples
    for msample in bamsexdata["M"]:
        if msample in clustersamples:
            cluster_fmpons["M"].append(msample)
    return cluster_fmpons


def write_output_file(outfileloc, ponsamples, ponoutdir, countsdir, clusnum):
    """Write a set of samples to use for the panel of normals.

    Parameters
    ----------
    outfileloc : str
        Path to write ouput file to
    ponsamples : list of str
        Samples to use for making the PanelOfNormals
    ponoutdir : str
        Path GATK4 should write PoN to
    countsdir : str
        Path to directory with counts file

    Returns
    -------
    file_written : bool
        True if output file has been succesfully written, False if not
    """
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write(f"#!/bin/bash\n")
            outfile.write(f"#SBATCH --job-name=solverd_pon_{clusnum}\n")
            outfile.write(f"#SBATCH --output=solverd_pon_{clusnum}.out\n")
            outfile.write(f"#SBATCH --error=solverd_pon_{clusnum}.err\n")
            outfile.write(f"#SBATCH --time=01:00:00\n")
            outfile.write(f"#SBATCH --cpus-per-task=1\n")
            outfile.write(f"#SBATCH --mem=7gb\n")
            outfile.write(f"#SBATCH --nodes=1\n")
            outfile.write(f"#SBATCH --open-mode=append\n")
            outfile.write(f"#SBATCH --export=NONE\n")
            outfile.write(f"#SBATCH --get-user-env=L\n\n")
            outfile.write("module load GATK/4.1.4.0-Java-8-LTS\n")
            outfile.write("module list\n\n")
            outfile.write("gatk CreateReadCountPanelOfNormals \\\n")
            for ponsample in ponsamples:
                outfile.write(f"\t-I {countsdir}{ponsample}.hdf5 \\\n")
            outfile.write(f"\t-O {ponoutdir}pon_{clusnum}.hdf5\n")
        file_written = True
    except IOError:
        print(f"Could not read output file {outfileloc}")
    finally:
        return file_written


def add_dir_slash(dirlocation):
    if not dirlocation.endswith("/"):
        return f"{dirlocation}/"
    return dirlocation


def main():
    fmcc_params = get_params()
    bam_sex_data = read_bamsex_data_2(fmcc_params["bam-to-sex"])
    s4_cluster_data = read_combineds4_data(fmcc_params["combined-s4"])

    outdir = add_dir_slash(fmcc_params["outdir"])
    ponoutdir = add_dir_slash(fmcc_params["pon-out"])
    countsdir = add_dir_slash(fmcc_params["read-counts-dir"])

    # Make female and male PaneOfNormals groups for each S4 cluster
    for clusternumber in s4_cluster_data:
        fmpon = make_fm_pon(s4_cluster_data[clusternumber], bam_sex_data)
        write_output_file(f"{outdir}fpon_{clusternumber}.sh", fmpon["F"], ponoutdir, countsdir, clusternumber)
        write_output_file(f"{outdir}mpon_{clusternumber}.sh", fmpon["M"], ponoutdir, countsdir, clusternumber)


if __name__ == "__main__":
    main()


# USAGE:
# python /groups/umcg-solve-rd/tmp01/umcg-mbeukers/solverd_scripts/fm_combined_clusters.py 
# -b /groups/umcg-solve-rd/tmp01/umcg-mbeukers/rd3_data/sex_to_bam.txt 
# -c /groups/umcg-solve-rd/tmp01/umcg-mbeukers/combined_s4_clusterwes/batch11.txt 
# -o /groups/umcg-solve-rd/tmp01/umcg-mbeukers/gatk4_freeze1/scripts/jobs/pon_jobs/batch11 
# -p /groups/umcg-solve-rd/tmp01/umcg-mbeukers/gatk4_freeze1/pon/batch11 
# -r /groups/umcg-solve-rd/tmp01/umcg-mbeukers/gatk4_freeze1/crc/batch11
