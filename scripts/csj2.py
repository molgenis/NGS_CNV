import os
import argparse


def add_dir_slash(dirtomodify):
    """Add a trailing '/' to a directory path

    Parameters
    ----------
    dirtomodify : str
        Path to directory

    Returns
    -------
    dirtomodify : str
        Path to directory with added trailing '/'
    """
    if not dirtomodify.endswith("/"):
        return f"{dirtomodify}/"
    return dirtomodify


def get_parameters():
    """Obtain and return set command line parameters

    Returns
    -------
    dict
        Command line parameter settings
    """
    # Make a list of valid generators
    valid_generators = ["ppi", "crc", "pon", "drc", "drc_n", "drc_p", "pdcr",
                        "cac", "ms", "ccrs", "pms"]

    sbatch_parameters = argparse.ArgumentParser()
    # Parameters related to sbatch job creation.
    sbatch_parameters.add_argument("-no", "--numofjobs", dest="numofjobs",
                                   required=True, type=int,
                                   help="Number of jobs to create.")
    sbatch_parameters.add_argument("-o", "--outdir", dest="outdir",
                                   required=True,
                                   help="Directory to write cluster jobs to.")

    # Parameters related to input directories.
    sbatch_parameters.add_argument("-ib", "--inaligments", dest="inalignments",
                                   help="Input directory with BAM/CRAM files.")
    sbatch_parameters.add_argument("-is", "--instandardized",
                                   dest="instandardized",
                                   help="Directory with standardized tsv "
                                        "count files.")
    sbatch_parameters.add_argument("-id", "--indenoised", dest="indenoised",
                                   help="Directory with denoised tsv count "
                                        "files.")
    sbatch_parameters.add_argument("-ia", "--inallelic", dest="inallelic",
                                   help="Directory with allelic tsv count "
                                        "files.")
    sbatch_parameters.add_argument("-ig", "--insegments", dest="insegments",
                                   help="Directory with segment files.")
    sbatch_parameters.add_argument("-il", "--intervallist",
                                   dest="intervallist",
                                   help="Interval list required for many "
                                        "processes by GATK4.")

    # Parameter related to the GATK4 commands themselves.
    sbatch_parameters.add_argument("-g", "--generate", dest="generate",
                                   required=True, choices=valid_generators,
                                   help="GATK4 tool sbatch scripts to "
                                        "generate.")
    sbatch_parameters.add_argument("-p", "--projectdir", dest="projectdir",
                                   required=True,
                                   help="Path to the project directory.")
    sbatch_parameters.add_argument("-i", "--indir", dest="indir",
                                   required=True,
                                   help="Directory with files to use for "
                                        "sbatch jobs.")
    sbatch_parameters.add_argument("-jo", "--joboutdir", dest="joboutdir",
                                   required=True,
                                   help="Output directory to use for each "
                                        "GATK command.")
    sbatch_parameters.add_argument("-r", "--genomeref", dest="genomeref",
                                   help="Genome reference in fasta format.")
    sbatch_parameters.add_argument("-rd", "--refdict", dest="refdict",
                                   help="Genome reference sequence dictionary "
                                        "file.")
    sbatch_parameters.add_argument("-pn", "--panelofnormals",
                                   dest="panelofnormals",
                                   help="GATK4 created Panel of Normals.")

    # Parameters related to the GATK4 commands themselves.
    sbatch_parameters.add_argument("-gm", "--gatkmem", dest="gatkmem",
                                   default="4g",
                                   help="Amount of memory for java VM running "
                                        "GATK4.")
    sbatch_parameters.add_argument("-bl", "--binlength", dest="binlength",
                                   default="0",
                                   help="Bin length to use for GATK commands "
                                        "where applicable.")
    sbatch_parameters.add_argument("-imr", "--intervalmergingrule",
                                   dest="intervalmergingrule",
                                   default="OVERLAPPING_ONLY",
                                   help="Interval merging rule to use for "
                                        "GATK jobs where applicable.")
    sbatch_parameters.add_argument("-mimp", "--minintervalmedianpercentile",
                                   dest="minimumintervalmedianpercentile",
                                   default="5.0",
                                   help="Minimum interval median percentile "
                                        "to use.")
    sbatch_parameters.add_argument("-mcl", "--minimumcontiglength",
                                   dest="minimumcontiglength",
                                   default="46709983",
                                   help="Minimum contig length to use for "
                                        "GATK4 commands where applicable.")
    sbatch_parameters.add_argument("-es", "--eigensamples", type=int,
                                   dest="eigensamples",
                                   help="Num of eigensamples for denoising.")

    # Parameters related to the SBATCH header.
    sbatch_parameters.add_argument("-j", "--jobname", dest="jobname",
                                   default="sbatchjob",
                                   help="Name to use for the sbatch jobs")
    sbatch_parameters.add_argument("-t", "--time", dest="time",
                                   default="04:00:00",
                                   help="Amount of time for each sbatch job")
    sbatch_parameters.add_argument("-c", "--cpus", dest="cpus", default="1",
                                   help="Number of cpus to use for each job")
    sbatch_parameters.add_argument("-m", "--mem", dest="mem", default="4gb",
                                   help="Amount of memory to use for each"
                                        "job (e.g. 4gb)")
    sbatch_parameters.add_argument("-n", "--nodes", dest="nodes", default="1",
                                   help="Number of nodes to use for each job")

    # Parameters related to loading of modules.
    sbatch_parameters.add_argument("-gv", "--gatkver", dest="gatkver",
                                   default="4.1.4.0-Java-8-LTS",
                                   help="GATK4 module to load and use")
    sbatch_parameters.add_argument("-pv", "--pyversion", dest="pythonversion",
                                   default="3.7.4-foss-2018b-v19.08.1",
                                   help="PythonPlus module to load and use")
    sbatch_parameters.add_argument("-rv", "-rversion", dest="rversion",
                                   default="3.6.1-foss-2018b-v19.11.1",
                                   help="RPlus module to load and use")

    # Other optional parameters for GATK4 parameters
    sbatch_parameters.add_argument("-op", "--optional", dest="optionalargs",
                                   default="",
                                   help="Optional arguments to set the for "
                                        "the GATK4 command to generate.")
    return vars(sbatch_parameters.parse_args())


def check_parameters(scriptparameters):
    """Check and return whether all provided parameters are set and ok.

    Parameters
    ----------
    scriptparameters : dict
        Script command line parameters

    Returns
    -------
    bool
        True if all parameters are ok, False if not
    """
    required_dir_params = ["projectdir", "indir", "outdir", "joboutdir"]
    required_gen_params = {"ppi": ["intervallist", "genomeref"],
                           "crc": ["intervallist"],
                           "drc": ["panelofnormals"],
                           "pdcr": ["refdict"],
                           "cac": ["intervallist", "genomeref"]}

    # Check whether the general parameters are ok
    dirparams_ok = check_dir_parameters(dict((x, scriptparameters[x])
                                             for x in required_dir_params
                                             if x in scriptparameters))

    # Check whether the required parameters for the specified generator are ok.
    gen_params_ok = True
    if scriptparameters["generate"] in required_gen_params:
        for genpar in required_gen_params[scriptparameters["generate"]]:
            if not os.path.isfile(scriptparameters[genpar]):
                gen_params_ok = False
                print(f"File {scriptparameters[genpar]} for parameter "
                      f"{genpar} does not seem to exist.")
    return dirparams_ok and gen_params_ok


def check_dir_parameters(dirparams):
    """Check and return whether provided directories via directory parameters
    exist.

    Parameters
    ----------
    dirparams : dict
        Directory parameters

    Returns
    -------
    dir_params_ok : bool
        True if all directories exist, False if not
    """
    dir_params_ok = True
    for dirpar in dirparams.values():
        if dirpar is not None:
            if dirpar != "":
                if not os.path.isdir(dirpar):
                    print(f"Directory {dirpar} does not exist")
                    dir_params_ok = False
            else:
                dir_params_ok = False
    return dir_params_ok


def divide_jobcommands_over_scripts(jobcommands, numofjobs):
    """Divide and return a list of job commands over a specified number of
    jobs.

    Parameters
    ----------
    jobcommands : list of str
        Job commands to divide
    numofjobs : int
        Number of sbatch jobs to divide commands over

    Returns
    -------
    split_coms : list of list of str
        Job commands divided over the requested number of job scripts
    """
    cycler = 0
    split_coms = [[] for x in range(numofjobs)]
    copy_list_jobcoms = list(jobcommands)
    while copy_list_jobcoms:
        cycle = cycler % numofjobs
        split_coms[cycle].append(copy_list_jobcoms.pop())
        cycler += 1
    return split_coms


def check_gatk_memory_to_job_memory(gatkmem, jobmem):
    """Check if the GATK4 java memory exceeds job memory. If so, scale down
    to job memory.

    Parameters
    ----------
    gatkmem : str
        Requested GATK4 java memory
    jobmem : str
        Requested job memory

    Returns
    -------
    gatk_mem : str
        Adjusted GATK4 java memory request
    """
    mem_measure_level = {"k": 0, "m": 1, "g": 2, "t": 3}
    gmem_num = int(gatkmem[0:-1])
    gmem_measure = gatkmem[-1].lower()
    jmem_num = int(jobmem[0:-2])
    jmem_measure = jobmem[-2].lower()

    # Check if the measure is the same
    if gmem_measure in mem_measure_level and jmem_measure in mem_measure_level:
        if mem_measure_level[gmem_measure] > mem_measure_level[jmem_measure]:
            gmem_measure = jmem_measure

            # Check if the gatk memory number is larger
            if gmem_num > jmem_num:
                gmem_num = jmem_num
    return f"{gmem_num}{gmem_measure}"


def generate_preprocess_intervals(intervallistloc, genomerefloc, imr,
                                  outdirloc, reqgatkmem):
    """Generate a command to pre-process an interval list

    Parameters
    ----------
    intervallistloc : str
        Path to the interval list to process
    genomerefloc : str
        Path to the genome reference fasta to use
    imr : str
        Interval merging rule to use for GATK4 commands
    outdirloc : str
        Path to the output directory for the GATK4 command
    reqgatkmem : str
        Requested GATK4 java VM memory

    Returns
    -------
    list of str
        GATK4 command to preprocess interval list
    """
    interval_prefix = get_file_prefix(intervallistloc.split("/")[-1])
    ppi_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
        f"PreprocessIntervals -L {intervallistloc} -R {genomerefloc} " \
        f"--bin-length 0 --interval-merging-rule {imr} " \
        f"-O {outdirloc}{interval_prefix}.processed.interval_list"
    return [ppi_com]


def generate_collect_read_counts(alnfiles, intervallistloc, imr, outdirloc,
                                 reqgatkmem):
    """Generate and return GATK4 commands to collect read counts from
    multiple BAM files.

    Parameters
    ----------
    alnfiles: list of str
        Alignment files to collect reads from
    intervallistloc: str
        Path to interval list to use for collecting read counts
    imr : str
        Interval merging rule to use in GATK4 commands
    outdirloc: str
        Path ot output directory for the GATK4 commands
    reqgatkmem : str
        Requested GATK4 java VM memory

    Returns
    -------
    crc_coms : list of str
        GATK4 CollectReadCounts commands
    """
    crc_coms = []
    for alnfile in alnfiles:
        hdf5_prefix = get_file_prefix(alnfile.split("/")[-1])
        crc_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
            f"CollectReadCounts -I {alnfile} -L {intervallistloc} " \
            f"--interval-merging-rule {imr} -O {outdirloc}{hdf5_prefix}.hdf5"
        crc_coms.append(crc_com)
    return crc_coms


def generate_panel_of_normals(hdf5files, mimp, outdirloc, reqgatkmem):
    """Generate and return a GATK4 command to create a panel of normals from
    HDF5 files.

    Parameters
    ----------
    hdf5files : list of str
        Path to HDF5 count files
    mimp : str
        Minimum interval median precentile to use in GATK4 commands
    outdirloc : str
        Path to output directory for GATK4 command
    reqgatkmem : str
        Requested GATK4 java VM memory

    Returns
    -------
    pon_coms : list of str
        GATK4 CreatePanelOfNormals command
    """
    pon_coms = []
    pon_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
        f"CreateReadCountPanelOfNormals \\\n"
    for hdf5file in hdf5files:
        pon_com += f"\t-I {hdf5file} \\\n"
    pon_com += f"\t--minimum-interval-median-percentile {mimp} \\\n"
    pon_com += f"\t-O {outdirloc}panel_of_normals.hdf5"
    pon_coms.append(pon_com)
    return pon_coms


def generate_denoise_read_counts(hdf5files, panelofnormals, outdirloc,
                                 reqgatkmem):
    """Generate and return GATK4 commands to denoise read counts against a
    panel of normals.

    Parameters
    ----------
    hdf5files : list of str
        HDF5 count files to denoise
    panelofnormals : str
        Panel of normals for the GATK commands
    outdirloc : str
        Path to output directory for GATK commands
    reqgatkmem : str
        Requested GATK4 java VM memory

    Returns
    -------
    drc_coms : list of str
        GATK 4 commands to denoise read counts
    """
    drc_coms = []
    for hdf5file in hdf5files:
        hdf5_prefix = get_file_prefix(hdf5file.split("/")[-1])
        drc_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
            f"DenoiseReadCounts -I {hdf5file} " \
            f"--count-panel-of-normals {panelofnormals} " \
            f"--standardized-copy-ratios {outdirloc}standardized/" \
            f"{hdf5_prefix}.standardized.tsv " \
            f"--denoised-copy-ratios {outdirloc}denoised/" \
            f"{hdf5_prefix}.denoised.tsv"
        drc_coms.append(drc_com)
    return drc_coms


def generate_denoise_read_counts_nopon(hdf5files, outdirloc, reqgatkmem):
    """Generate and return GATK4 commands to denoise read counts without a
    panel of normals.

    Parameters
    ----------
    hdf5files : list of str
        HDF5 count files to denoise.
    outdirloc : str
        Path to output directory for GATK4 commands.
    reqgatkmem : str
        Requested GATK4 java VM memory.

    Returns
    -------
    drc_coms : list of str
        GATK4 command to denoise read counts without a PoN.
    """
    drc_coms = []
    for hdf5file in hdf5files:
        hdf5_prefix = get_file_prefix(hdf5file.split("/")[-1])
        drc_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
            f"DenoiseReadCounts -I {hdf5file} " \
            f"--standardized-copy-ratios {outdirloc}standardized/" \
            f"{hdf5_prefix}.standardized.tsv " \
            f"--denoised-copy-ratios {outdirloc}denoised/" \
            f"{hdf5_prefix}.denoised.tsv"
        drc_coms.append(drc_com)
    return drc_coms


def generate_denoise_read_counts_pcs(hdf5files, panelofnormals, numofpcs,
                                     outdirloc, reqgatkmem):
    """Generate and return GATK4 commands ot denoise read counts with a panel
    of normals and a set number of principal components / eigensamples.

    Parameters
    ----------
    hdf5files : list of str
        HDF5 count files to denoise.
    panelofnormals : str
        Panel of normals to use.
    numofpcs : int
        Number of PCs/eigensamples to use for denoising.
    outdirloc : str
        Path to output directory for GATK4 commands.
    reqgatkmem : str
        Requested GATK4 java VM memory.

    Returns
    -------
    drc_coms : list of str
        GATK4 commands to denoise read counts with a set number of PCs.
    """
    drc_coms = []
    for hdf5file in hdf5files:
        hdf5_prefix = get_file_prefix(hdf5file.split("/")[-1])
        drc_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
            f"DenoiseReadCounts -I {hdf5file} " \
            f"--count-panel-of-normals {panelofnormals} " \
            f"--standardized-copy-ratios {outdirloc}standardized/" \
            f"{hdf5_prefix}.standardized.tsv " \
            f"--denoised-copy-ratios {outdirloc}denoised/" \
            f"{hdf5_prefix}.denoised.tsv " \
            f"--number-of-eigensamples {numofpcs}"
        drc_coms.append(drc_com)
    return drc_coms


def generate_plot_denoised_copy_ratios(tsv_standardized, tsv_denoised,
                                       sequencedict, mcl, outdirloc,
                                       reqgatkmem):
    """Generate and return GATK4 commands for plotting denoised copy ratios.

    Parameters
    ----------
    tsv_standardized : list of str
        List of standardized hdf5 count files
    tsv_denoised : list of str
        List of denoised hdf5 count files
    sequencedict : str
        Genome reference sequence dictionary
    mcl : str
        Minimum contig length for GATK4 commands
    outdirloc : str
        Path to output directory for GATK4 command
    reqgatkmem : str
        Request GATK4 java VM memory

    Returns
    -------
    pdcr_coms : list of str
        GATK4 plot denoised copy ratio commands
    """
    # Sort the tsv_standardized and tsv_denoised lists to ensure the same order
    tsv_standardized.sort()
    tsv_denoised.sort()

    pdcr_coms = []
    for x in range(0, len(tsv_standardized)):
        tsv_prefix = get_file_prefix(tsv_standardized[x].split("/")[-1])
        pdcr_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
            f"PlotDenoisedCopyRatios " \
            f"--standardized-copy-ratios {tsv_standardized[x]} " \
            f"--denoised-copy-ratios {tsv_denoised[x]} " \
            f"--sequence-dictionary {sequencedict} " \
            f"--minimum-contig-length {mcl} --output {outdirloc} " \
            f"--output-prefix {tsv_prefix}"
        pdcr_coms.append(pdcr_com)
    return pdcr_coms


def generate_collect_allelic_counts(alnfiles, intervallistloc, genomereference,
                                    outdirloc, reqgatkmem):
    """Generate and return GATK4 commands for performing CollectAllelicCounts.

    Parameters
    ----------
    alnfiles : list of str
        List of alignment files (BAM/CRAM) to collect allelic counts from
    intervallistloc : str
        Path to intervallist to use
    genomereference : str
        Path to genome reference to use
    outdirloc : str
        Path to output directory to write
    reqgatkmem : str
        Requested GATK4 java VM memory

    Returns
    -------
    cac_coms : list of str
        GATK4 CollectAllelicCounts commands for supplied alignment files
    """
    cac_coms = []
    for alnfile in alnfiles:
        tsv_prefix = get_file_prefix(alnfile.split("/")[-1])
        cac_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
            f"CollectAllelicCounts -L {intervallistloc} -I {alnfile} " \
            f"-R {genomereference} -O {outdirloc}{tsv_prefix}.allelecounts.tsv"
        cac_coms.append(cac_com)
    return cac_coms


def generate_model_segments_s(tsv_denoised, tsv_allelic, outdir, reqgatkmem):
    """Generate and return GATK4 commands for performing ModelSegments.

    This method generates GATK4 ModelSegments only using sample denoised and
    allelic count files as input.

    Parameters
    ----------
    tsv_denoised : list of str
        Denoised sample tsv count files.
    tsv_allelic : list of str
        Allelic sample tsv count files.
    outdir : str
        Directory to write output to.
    reqgatkmem : str
        Requested GATK4 java VM memory.

    Returns
    -------
    ms_coms : list of str
        GATK4 ModelSegments commands.
    """
    ms_coms = []
    for x in range(0, len(tsv_denoised)):
        ms_prefix = get_file_prefix(tsv_denoised[x].split("/")[-1])
        ms_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
            f"ModelSegments --denoised-copy-ratios {tsv_denoised[x]} " \
            f"--allelic-counts {tsv_allelic[x]} --output {outdir} " \
            f"--output-prefix {ms_prefix}"
        ms_coms.append(ms_com)
    return ms_coms


def generate_model_segments(tsv_denoised, tsv_samplealleliccounts,
                            tsv_refalleliccounts, outdir, reqgatkmem):
    """Generate and return GATK4 commands for performing ModelSegments.

    Parameters
    ----------
    tsv_denoised : list of str
        Denoised sample counts files
    tsv_samplealleliccounts : list of str
        Sample allelic count files
    tsv_refalleliccounts : list of str
        Reference allelic count files
    outdir : str
        Path to GATK4 job output directory
    reqgatkmem : str
        Requested GATK4 java VM memory

    Returns
    -------
    ms_coms : list of str
        GATK4 ModelSegments commands
    """
    ms_coms = []
    for x in range(0, len(tsv_denoised)):
        ms_prefix = get_file_prefix(tsv_denoised[x].split("/")[-1])
        ms_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
            f"ModelSegments --denoised-copy-ratios {tsv_denoised[x]} " \
            f"--allelic-counts {tsv_samplealleliccounts[x]} " \
            f"--normal-allelic-counts {tsv_refalleliccounts[x]} " \
            f"--output {outdir} --output-prefix {ms_prefix}"
        ms_coms.append(ms_com)
    return ms_coms


def generate_call_copy_ratio_segments(segfiles, outdirloc, reqgatkmem):
    """Generate and return GATK4 commands for performing CallCopyRatioSegments.

    Parameters
    ----------
    segfiles : list of str
        Paths to segment files
    outdirloc : str
        Path to output directory for GATK4 commands
    reqgatkmem : str
        Requested GATK4 java VM memory

    Returns
    -------
    ccrs_coms : list of str
        GATK4 CallCopyRatioSegments commands
    """
    ccrs_coms = []
    for segfile in segfiles:
        outname_prefix = get_file_prefix(segfile.split("/")[-1])
        ccrs_com = f"gatk --java-options \"-Xmx{reqgatkmem}\" " \
            f"CallCopyRatioSegments --input {segfile} " \
            f"--output {outdirloc}{outname_prefix}.called.seg"
        ccrs_coms.append(ccrs_com)
    return ccrs_coms


def generate_plot_modeled_segments(tsv_denoised, tsv_allelic, segfiles,
                                   sequencedict, mcl, outdir, reqgatkmem):
    """Generate and return GATK4 commands for performing PlotModeledSegments.

    Parameters
    ----------
    tsv_denoised : list of str
        Denoised counts files
    tsv_allelic : list of str
        Allelic counts files
    segfiles : list of str
        Segment files
    sequencedict : str
        Genome reference sequence dictionary
    mcl : str
        Minimum contig length for GATK4 commands
    outdir : str
        Path to output directory for GATK4 commands
    reqgatkmem : str
        Requested GATK4 java VM memory

    Returns
    -------
    pms_coms : list of str
        GATK4 PlotModelledSegments commands
    """
    pms_coms = []
    for x in range(0, len(tsv_denoised)):
        outname_prefix = get_file_prefix(tsv_denoised[x].split("/")[-1])
        pms_com = f"gatk --java-options\"-Xmx{reqgatkmem}\" " \
            f"PlotModeledSegments --denoised-copy-ratios {tsv_denoised[x]} " \
            f"--allelic-counts {tsv_allelic[x]} --segments {segfiles[x]} " \
            f"--sequence-dictionary {sequencedict} " \
            f"--minimum-contig-length {mcl}  --output {outdir} " \
            f"--output-prefix {outname_prefix}"
    return pms_coms


def get_required_files(filesdir, fileext):
    """Extract and return the list of HDF5 count files in the provided
    directory.

    Parameters
    ----------
    filesdir : str
        Directory containing HDF5 count files within the project directory
    fileext : str or tuple of str
        File extension(s) to select specific files

    Returns
    -------
    directory_files : list of str
        Paths to input files
    """
    directory_files = []
    if os.path.isdir(filesdir):
        directory_files = os.listdir(filesdir)
        directory_files = [f"{filesdir}{x}" for x in directory_files
                           if x.endswith(fileext)]
    return directory_files


def get_file_prefix(filename):
    """Obtain and return the prefix of a filename.

    Parameters
    ----------
    filename : str
        Complete file name to extract prefix from

    Returns
    -------
    str
        Prefix of supplied file
    """
    return filename.split(".")[0]


def generate_sbatch_header(sbatch_values):
    """Generate and return the SBATCH header.

    Parameters
    ----------
    sbatch_values : dict
        Settings for the SBATCH job header

    Returns
    -------
    sbatchheader : str
        Header with the set parameters
    """
    jobname = sbatch_values["jobname"]
    jobtime = sbatch_values["time"]
    jobcpus = sbatch_values["cpus"]
    jobmem = sbatch_values["mem"]
    jobnodes = sbatch_values["nodes"]

    sbatchheader = f"#!/bin/bash\n" \
        f"#SBATCH --job-name={jobname}\n" \
        f"#SBATCH --output={jobname}.out\n" \
        f"#SBATCH --error={jobname}.err\n" \
        f"#SBATCH --time={jobtime}\n" \
        f"#SBATCH --cpus-per-task={jobcpus}\n" \
        f"#SBATCH --mem={jobmem}\n" \
        f"#SBATCH --nodes={jobnodes}\n" \
        f"#SBATCH --open-mode=append\n" \
        f"#SBATCH --export=NONE\n" \
        f"#SBATCH --get-user-env=L\n\n"
    return sbatchheader


def make_sbatch_jobs(jobname, sbatch_settings, num_of_jobs, jobcommands,
                     gatkver, rversion, outdir):
    """Make and write a specified number of sbatch jobs to separate sbatch
    job files.

    Parameters
    ----------
    jobname : str
        Name suffix for sbatch jobs
    sbatch_settings : dict
        Settings to place in sbatch header
    num_of_jobs : int
        Number of sbatch jobs to create
    jobcommands : list of str
        Job commands to divide and place in sbatch job scripts
    gatkver : str
        GATK4 version to use
    rversion : str
        RPlus module version to load and use
    outdir : str
        Directory to write the sbatch jobs to
    """
    # Check whether the number of commands is smaller than the number of
    # requested jobs.
    num_of_jobcoms = len(jobcommands)
    if num_of_jobcoms < num_of_jobs:
        print(f"The number of job commands ({num_of_jobcoms}) "
              f"is smaller than the number of requested "
              f"sbatch scripts ({num_of_jobs}).")
        print(f"The number of sbatch scripts will be {num_of_jobcoms}")
        num_of_jobs = num_of_jobcoms

    # Divide the sbatch job command over the number of sbatch jobs to create
    divided_jobcoms = divide_jobcommands_over_scripts(jobcommands, num_of_jobs)

    # Make the sbatch job scripts
    for x in range(1, num_of_jobs + 1):
        outpath = f"{outdir}{jobname}_{x}.sh"
        sbatch_settings["jobname"] = f"{jobname}_{x}"
        sbatch_header = generate_sbatch_header(sbatch_settings)
        create_sbatch_script(sbatch_header, divided_jobcoms[x-1], gatkver,
                             rversion, outpath)


def create_sbatch_script(header, commands, gatkver, rversion, outfilepath):
    """Write a single sbatch script file.

    Parameters
    ----------
    header : str
        Sbatch header to add
    commands : list of str
        Command to place in sbatch script
    gatkver : str
        Version of GATK4 module to use
    rversion : str
        Version of RPlus module to use
    outfilepath : str
        Output path to write sbatch script file to
    """
    try:
        with open(outfilepath, "w") as scriptfile:
            scriptfile.write(header)
            scriptfile.write(f"module load GATK/{gatkver}\n")
            scriptfile.write(f"module load RPlus/{rversion}\n")
            scriptfile.write("module list\n\n")
            scriptfile.write("\n".join(commands))
    except IOError:
        print(f"Could not write script file {outfilepath} :\'(")


def add_optional_parameters(jobcommands, optparams):
    """Adds specified optional parameters to created job commands.

    Parameters
    ----------
    jobcommands : list of str
        Job commands to add optional parameters to.
    optparams : str
        Optional parameters to add to job commands.

    Returns
    -------
    list of str
        Job commands to write to job script files
    """
    return [jobcom + f" {optparams}" for jobcom in jobcommands]


if __name__ == "__main__":
    sbatch_gen_params = get_parameters()

    projectdir = add_dir_slash(sbatch_gen_params["projectdir"])
    inputdir = add_dir_slash(sbatch_gen_params["indir"])
    sbatch_outdir = add_dir_slash(sbatch_gen_params["outdir"])
    gatkjob_outdir = add_dir_slash(sbatch_gen_params["joboutdir"])
    requested_gatkmem = check_gatk_memory_to_job_memory(sbatch_gen_params["gatkmem"], sbatch_gen_params["mem"])
    optparams = sbatch_gen_params["optionalargs"]
    job_commands = []

    # Check whether the parameters are ok and decide which generator to run
    if check_parameters(sbatch_gen_params):
        indirfiles = []

        # Generate sbatch command for PreprocessIntervals
        if sbatch_gen_params["generate"] == "ppi":
            job_commands = generate_preprocess_intervals(sbatch_gen_params["intervallist"],
                                                         sbatch_gen_params["genomeref"],
                                                         sbatch_gen_params["intervalmergingrule"], gatkjob_outdir,
                                                         requested_gatkmem)

        # Generate sbatch commands for CollectReadCounts
        elif sbatch_gen_params["generate"] == "crc":
            indirfiles = get_required_files(inputdir, (".bam", ".cram"))
            job_commands = generate_collect_read_counts(indirfiles, sbatch_gen_params["intervallist"],
                                                        sbatch_gen_params["intervalmergingrule"], gatkjob_outdir,
                                                        requested_gatkmem)

        # Generate sbatch commands for CreatePanelOfNormals
        elif sbatch_gen_params["generate"] == "pon":
            indirfiles = get_required_files(inputdir, ".hdf5")
            job_commands = generate_panel_of_normals(indirfiles, sbatch_gen_params["minimumintervalmedianpercentile"],
                                                     gatkjob_outdir, requested_gatkmem)

        # Generate sbatch commands for DenoiseReadCounts
        elif sbatch_gen_params["generate"] == "drc":
            indirfiles = get_required_files(inputdir, ".hdf5")
            job_commands = generate_denoise_read_counts(indirfiles, sbatch_gen_params["panelofnormals"], gatkjob_outdir,
                                                        requested_gatkmem)

        # Generate sbatch commands for DenoiseReadCounts without a PoN.
        elif sbatch_gen_params["generate"] == "drc_n":
            indirfiles = get_required_files(inputdir, ".hdf5")
            job_commands = generate_denoise_read_counts_nopon(indirfiles,
                                                              gatkjob_outdir,
                                                              requested_gatkmem)

        # Generate sbatch commands for DenoiseReadCounts with eigensamples
        elif sbatch_gen_params["generate"] == "drc_p":
            indirfiles = get_required_files(inputdir, ".hdf5")
            job_commands = generate_denoise_read_counts_pcs(indirfiles, sbatch_gen_params["panelofnormals"],
                                                            sbatch_gen_params["eigensamples"], gatkjob_outdir,
                                                            requested_gatkmem)

        # Generate sbatch commands for PlotDenoisedReadCounts
        elif sbatch_gen_params["generate"] == "pdcr":
            standardized_infiles = get_required_files(f"{inputdir}standardized/", ".tsv")
            denoised_infiles = get_required_files(f"{inputdir}denoised/", ".tsv")
            job_commands = generate_plot_denoised_copy_ratios(standardized_infiles, denoised_infiles,
                                                              sbatch_gen_params["refdict"],
                                                              sbatch_gen_params["minimumcontiglength"], gatkjob_outdir,
                                                              requested_gatkmem)

        # Generate sbatch commands for CollectAllelicCounts
        elif sbatch_gen_params["generate"] == "cac":
            indirfiles = get_required_files(inputdir, (".bam", ".cram"))
            job_commands = generate_collect_allelic_counts(indirfiles, sbatch_gen_params["intervallist"],
                                                           sbatch_gen_params["genomeref"], gatkjob_outdir,
                                                           requested_gatkmem)

        # Generate sbatch command for ModelSegments
        elif sbatch_gen_params["generate"] == "ms":
            denoised_dir = add_dir_slash(sbatch_gen_params["indenoised"])
            denoised_infiles = get_required_files(denoised_dir, ".tsv")
            allelic_dir = add_dir_slash(sbatch_gen_params["inallelic"])
            allelic_infiles = get_required_files(allelic_dir, ".tsv")

            denoised_infiles.sort()
            allelic_infiles.sort()

            job_commands = generate_model_segments_s(denoised_infiles,
                                                     allelic_infiles,
                                                     gatkjob_outdir,
                                                     requested_gatkmem)

        # Generate sbatch commands for CallCopyRatioSegments
        elif sbatch_gen_params["generate"] == "ccrs":
            indirfiles = get_required_files(inputdir, ".cr.seg")
            job_commands = generate_call_copy_ratio_segments(indirfiles,
                                                             gatkjob_outdir,
                                                             requested_gatkmem)

        # Generate sbatch commands for PlotModeledSegments
        elif sbatch_gen_params["generate"] == "pms":
            # Fix the directories by adding a '/'.
            denoised_dir = add_dir_slash(sbatch_gen_params["indenoised"])
            allelic_dir = add_dir_slash(sbatch_gen_params["inallelic"])
            segments_dir = add_dir_slash(sbatch_gen_params["insegments"])

            # Get the required files.
            denoised_infiles = get_required_files(denoised_dir, ".tsv")
            allelic_infiles = get_required_files(allelic_dir, ".tsv")
            segment_infiles = get_required_files(segments_dir, ".cr.seg")
            refdict = sbatch_gen_params["refdict"]
            mcl = sbatch_gen_params["minimumcontiglength"]
            gatkmem = sbatch_gen_params["gatkmem"]

            # Sort the input files they are in the same order.
            denoised_infiles.sort()
            allelic_infiles.sort()
            segment_infiles.sort()

            job_commands = generate_plot_modeled_segments(denoised_infiles,
                                                          allelic_infiles,
                                                          segment_infiles,
                                                          refdict, mcl,
                                                          gatkjob_outdir,
                                                          gatkmem)

    # Check whether to add any optional parameters
    if optparams is not None and optparams != "":
        job_commands = add_optional_parameters(job_commands, optparams)

    # Generate the sbatch scripts with the constructed commands
    make_sbatch_jobs(sbatch_gen_params["jobname"], sbatch_gen_params,
                     sbatch_gen_params["numofjobs"], job_commands,
                     sbatch_gen_params["gatkver"],
                     sbatch_gen_params["rversion"], sbatch_outdir)
