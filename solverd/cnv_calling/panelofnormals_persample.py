import os
import argparse


def get_params():
    pon_args = argparse.ArgumentParser()
    pon_args.add_argument("-i", "--infile", required=True, help="Path to Lennart file with samples and cluster numbers")
    pon_args.add_argument("-b", "--batchtokit", required=True, dest="batchtokit", help="")
    pon_args.add_argument("-c", "--countsdir", required=True, dest="countsdir", help="Path to directory containg counts files from CollectReadCounts")
    pon_args.add_argument("-k", "--kitname", required=True, dest="kitname", help="Name of the kit")
    pon_args.add_argument("-o", "--outdir", required=True, dest="outdir", help="Path to write PanelOfNormals jobs to")
    pon_args.add_argument("-j", "--jobout", required=True, dest="jobout", help="Path for the job to write PoN to")
    pon_args.add_argument("--persample", dest="persample", action="store_true")
    return vars(pon_args.parse_args())


def add_dir_slash(dirloc):
    if not dirloc.endswith("/"):
        return f"{dirloc}/"
    return dirloc


def get_count_files(countsdirloc):
    """Return paths to counts files per sample name.

    Returns
    -------
    count_files : dict
        Paths to HDF5 count files per sample name
    """
    count_files = {}
    countsdirloc = add_dir_slash(countsdirloc)
    dirfiles = os.listdir(countsdirloc)
    hdf5_files = [hdffile for hdffile in dirfiles if hdffile.endswith(".hdf5")]
    for hdffile in hdf5_files:
        count_files[hdffile.split(".")[0]] = f"{countsdirloc}{hdffile}"
    return count_files


def get_batch_to_kit(batchtokitloc):
    batchtokit = {}
    try:
        with open(batchtokitloc, 'r')as btkfile:
            for btkline in btkfile:
                btkdata = btkline.strip().split("\t")
                batchtokit[btkdata[1]] = btkdata[0]
    except IOError:
        print(f"Could not open {batchtokitloc}")
    finally:
        return batchtokit


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


def determine_pons(clusterfileloc):
    """Read the Lennart cluster file.

    Parameters
    ----------
    clusterfileloc : str
        Path to Lennart sample cluster file

    Returns
    -------
    """
    batch_clusters = {}
    try:
        with open(clusterfileloc, 'r') as clusterfile:
            for clusterline in clusterfile:
                if clusterline.startswith("["):
                    batch_clusters = get_cluster_samples(clusterline.strip()[1:-1], batch_clusters)
    except IOError:
        print(f"Could not read Lennart cluster file {clusterfileloc}")
    finally:
        return batch_clusters


def make_single_pon(batchname, batchclusters, countfiles, joboutdirloc, outdirloc):
    outdirloc = add_dir_slash(outdirloc)
    joboutdirloc = add_dir_slash(joboutdirloc)
    for batchnum in batchclusters:
        outfilename = f"{outdirloc}createpon_{batchname}_{batchnum}.sh"
        try:
            with open(outfilename, 'w') as outfile:
                outfile.write("#!/bin/bash\n")
                outfile.write(f"#SBATCH --job-name=pon_{batchname}_{batchnum}\n")
                outfile.write(f"#SBATCH --output=pon_{batchname}_{batchnum}.out\n")
                outfile.write(f"#SBATCH --error=pon_{batchname}_{batchnum}.err\n")
                outfile.write("#SBATCH --time=02:00:00\n")
                outfile.write("#SBATCH --cpus-per-task=1\n")
                outfile.write("#SBATCH --mem=4gb\n")
                outfile.write("#SBATCH --nodes=1\n")
                outfile.write("#SBATCH --open-mode=append\n")
                outfile.write("#SBATCH --export=NONE\n")
                outfile.write("#SBATCH --get-user-env=L\n\n")

                outfile.write("module load GATK/4.1.4.1-Java-8-LTS\n")
                outfile.write("module load RPlus/3.6.1-foss-2018b-v19.11.1\n")
                outfile.write("module list\n\n")

                outfile.write("gatk CreateReadCountPanelOfNormals \\\n")
                for samplename in batchclusters[batchnum]:
                    if samplename in countfiles:
                        outfile.write(f"\t-I {countfiles[samplename]} \\\n")
                outfile.write("\t--minimum-interval-median-percentile 5.0 \\\n")
                outfile.write(f"\t-O {joboutdirloc}pon_{batchname}_{batchnum}.hdf5\n")
        except:
            print("AAP")


def main():
    pon_params = get_params()
    countfiles = get_count_files(pon_params["countsdir"])
    batch_to_kit = get_batch_to_kit(pon_params["batchtokit"])
    batchname = batch_to_kit[pon_params["kitname"]]
    batch_clusters = determine_pons(pon_params["infile"])
    make_single_pon(batchname, batch_clusters, countfiles, pon_params["jobout"], pon_params["outdir"])


if __name__ == "__main__":
    main()
