import argparse


def get_params():
    udpon_args = argparse.ArgumentParser()
    udpon_args.add_argument("-c", "--s4-cluster", type=str, required=True, dest="s4-cluster", help="Path to S4 ClusterWES file")
    udpon_args.add_argument("-n", "--name", type=str, required=True, dest="name", help="Prefix name to use (like batch11)")
    udpon_args.add_argument("-o", "--outdir", type=str, required=True, dest="outdir", help="Path to write output files to")
    udpon_args.add_argument("-r", "--read-counts-dir", type=str, required=True, dest="read-counts-dir", help="Path to directory with read count files")
    udpon_args.add_argument("-u", "--ud-detsex-file", type=str, required=True, dest="ud-detsex-file", help="Path to detsex file")
    return vars(udpon_args.parse_args())


def read_ud_detsex(uddetsexfileloc):
    ud_detsex_data = {"F": [], "M": []}
    try:
        with open(uddetsexfileloc, 'r') as detsexfile:
            for detsexline in detsexfile:
                detsexdata = detsexline.strip().split(" ")
                if detsexdata[1] == "F":
                    ud_detsex_data["F"].append(detsexdata[0][:-1])
                if detsexdata[1] == "M" or detsexdata[1] == "UD":
                    ud_detsex_data["M"].append(detsexdata[0][:-1])
    except IOError:
        print(f"Could not read {uddetsexfileloc}")
    finally:
        return ud_detsex_data


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


def write_udsample_to_pon(outfilepath, s4clustersamples, readcountdir, ponpath, udsamples):
    file_written = False
    try:
        with open(outfilepath, 'w') as outfile:
            for s4sample in s4clustersamples:
                if s4sample in udsamples:
                    outfile.write(f"{readcountdir}{s4sample}.hdf5\t{ponpath}\n")
        file_written = True
    except IOError:
        print(f"Could not write output file {outfilepath}")
    finally:
        return file_written


def main():
    udpon_params = get_params()
    s4data = read_combineds4_data(udpon_params["s4-cluster"])
    uddetsexdata = read_ud_detsex(udpon_params["ud-detsex-file"])

    readcountdir = udpon_params["read-counts-dir"]
    readcountdir = f"{readcountdir}/" if not readcountdir.endswith("/") else readcountdir
    outdirloc = udpon_params["outdir"]
    outdirloc = f"{outdirloc}/" if not outdirloc.endswith("/") else outdirloc
    prefixname = udpon_params["name"]

    for clusternum in s4data:
        f_udpon_outpath = f"{outdirloc}{prefixname}_{clusternum}_udfsamples.txt"
        m_udpon_outpath = f"{outdirloc}{prefixname}_{clusternum}_udmsamples.txt"
        f_ponfilepath = f"/groups/umcg-solve-rd/tmp01/umcg-mbeukers/gatk4_freeze1/pon/{prefixname}/{prefixname}_fpon_{clusternum}.hdf5"
        m_ponfilepath = f"/groups/umcg-solve-rd/tmp01/umcg-mbeukers/gatk4_freeze1/pon/{prefixname}/{prefixname}_mpon_{clusternum}.hdf5"
        wrote_f = write_udsample_to_pon(f_udpon_outpath, s4data[clusternum], readcountdir, f_ponfilepath, uddetsexdata["F"])
        wrote_m = write_udsample_to_pon(m_udpon_outpath, s4data[clusternum], readcountdir, m_ponfilepath, uddetsexdata["M"])
        print(f"Wrote male UD samples to PoN output file?: {wrote_m}")
        print(f"Wrote female UD samples to PoN output file?: {wrote_f}")


if __name__ == "__main__":
    main()
