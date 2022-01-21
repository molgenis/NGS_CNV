import os
import argparse


def get_params():
    """Define, receive and return CLI parameter values."""
    outchoices = ["sample", "summary", "all"]
    noct_args = argparse.ArgumentParser()
    noct_args.add_argument("-c", "--calls-dir", type=str, dest="calls-dir", help="Path to directory containing the CNV calls")
    noct_args.add_argument("-o", "--output-type", type=str, choices=outchoices, dest="output-type", help="Type of output to show")
    noct_args.add_argument("-od", "--out-dir", type=str, dest="out-dir", help="")
    noct_args.add_argument("-f", "--files-only", dest="files-only", action="store_true", help="")
    noct_args.add_argument("-p", "--prefix", type=str, dest="prefix", help="")
    return vars(noct_args.parse_args())


def get_num_of_call_types(callsfileloc):
    """Get and return the call types for a single sample."""
    call_type_counts = {'+':0, '0':0, '-':0}
    try:
        with open(callsfileloc, 'r') as callsfile:
            next(callsfile)
            for fileline in callsfile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[5] in call_type_counts.keys():
                    call_type_counts[filelinedata[5]] += 1
    except IOError:
        print("")
    finally:
        return call_type_counts


def summarize_samples(sample_counts):
    """Summarize and return calls over all samples."""
    summarized_counts = {'+':0, '0':0, '-':0}
    for samplename in sample_counts:
        summarized_counts['+'] += sample_counts[samplename]['+']
        summarized_counts['0'] += sample_counts[samplename]['0']
        summarized_counts['-'] += sample_counts[samplename]['-']
    return summarized_counts


def summarize_call_types(callsdir, callsfiles):
    summarized_counts = {'+':0, '0':0, '-':0}
    for ccrsfile in callsfiles:
        samplecounts = get_num_of_call_types(f"{callsdir}{ccrsfile}")
        summarized_counts['+'] += samplecounts['+']
        summarized_counts['0'] += samplecounts['0']
        summarized_counts['-'] += samplecounts['-']
    return summarized_counts


def perform_all_samples(callsdir, calls_files, fileoutpath, fileonly):
    sample_numbers = {}
    for callsample in calls_files:
        samplename = callsample.split(".")[0]
        sample_ctnums = get_num_of_call_types(f"{callsdir}{callsample}")
        sample_numbers[samplename] = sample_ctnums

        if not fileonly:
            # Display the sample numbers
            print(f"Sample: {samplename}")
            print("Number of gain calls: " +str(sample_ctnums['+']))
            print("Number of neutrals: " +str(sample_ctnums['0']))
            print("Number of loss calls: " +str(sample_ctnums['-']))
    write_all_samples(fileoutpath, sample_numbers)
    return sample_numbers


def write_all_samples(outfilepath, samplenums):
    try:
        with open(outfilepath, 'w') as outfile:
            outfile.write("Samples\tDuplication\tDeletion\n")
            for samplename in samplenums:
                outfile.write(samplename+ "\t" +str(samplenums[samplename]['+'])+ "\t" +str(samplenums[samplename]['-'])+ "\n")
    except IOError:
        print("Could not write output file")


def perform_all_summary(callsdir, callsfiles, samplenumbers, fileoutpath, fileonly):
    if len(samplenumbers) > 0:
        summarized_counts = summarize_samples(samplenumbers)
        if not fileonly:
            print("")
            print("")
    else:
        summarized_counts = summarize_call_types(callsdir, callsfiles)

    if not filelony:
        # Display the summary numbers
        print("[-CCRS Summary-]")
        print(f"Number of samples: {len(calls_files)}")
        print("Number of gain calls: " +str(summarized_counts['+']))
        print("Number of neutrals: " +str(summarized_counts['0']))
        print("Number of loss calls: " +str(summarized_counts['-']))
        print("Average gain calls per sample: " +str(round(summarized_counts['+']/len(calls_files), 3)))
        print("Average neutrals per sample: " +str(round(summarized_counts['0']/len(calls_files), 3)))
        print("Average loss calls per sample: " +str(round(summarized_counts['-']/len(calls_files), 3)))
    write_all_summary(fileoutpath, callsfiles, summarized_counts)


def write_all_summary(outfilepath, calls_files, summarized_counts):
    try:
        with open(outfilepath, 'w') as outfile:
            outfile.write("[-CCRS Summary-]\n")
            outfile.write(f"Number of samples: {len(calls_files)}\n")
            outfile.write("Number of gain calls: " +str(summarized_counts['+'])+ "\n")
            outfile.write("Number of neutrals: " +str(summarized_counts['0'])+ "\n")
            outfile.write("Number of loss calls: " +str(summarized_counts['-'])+ "\n")
            outfile.write("Average gain calls per sample: " +str(round(summarized_counts['+']/len(calls_files), 3))+ "\n")
            outfile.write("Average neutrals per sample: " +str(round(summarized_counts['0']/len(calls_files), 3))+ "\n")
            outfile.write("Average loss calls per sample: " +str(round(summarized_counts['-']/len(calls_files), 3))+ "\n")
    except IOError:
        print("")


def main():
    """Do the main work."""
    noct_params = get_params()

    callsdir = noct_params["calls-dir"]
    callsdir = f"{callsdir}/" if not callsdir.endswith("/") else callsdir
    calls_files = [x for x in os.listdir(callsdir) if x.endswith(".igv.seg")]

    outdir = noct_params["out-dir"]
    outdir = f"{outdir}/" if not outdir.endswith("/") else outdir
    outprefix = noct_params["prefix"]

    sample_numbers = {}
    summary_numbers = {}

    # Display number of call types per sample
    if noct_params["output-type"] == "sample" or noct_params["output-type"] == "all":
        outfileloc = f"{outdir}{outprefix}_dupdel.txt"
        sample_numbers = perform_all_samples(callsdir, calls_files, outfileloc, noct_params["files-only"])

    # Display a summary of call types over all samples
    if noct_params["output-type"] == "summary" or noct_params["output-type"] == "all":
        outfileloc = f"{outdir}{outprefix}_dupdel_summary.txt"
        perform_all_summary(callsdir, calls_files, sample_numbers, outfileloc, noct_params["files-only"])


if __name__ == "__main__":
    main()
