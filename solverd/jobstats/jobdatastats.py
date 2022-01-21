import argparse
import statistics


def get_params():
    jds_args = argparse.ArgumentParser()
    jds_args.add_argument("-i", "--infile", dest="infile", required=True, help="Path to the combined job out stats datafile.")
    jds_args.add_argument("-j", "--jobname", dest="jobname", help="Name of the job that generated the statistics")
    return vars(jds_args.parse_args())


def read_stats(statsfile):
    statsdata = {"JobID": [],
                 "Elapsed": [],
                 "AllocCPUS": [],
                 "AveCPU": [],
                 "ReqMem": [],
                 "MaxVMSize": [],
                 "MaxRSS": [],
                 "MaxDiskRead": [],
                 "MaxDiskWrite": []}
    try:
        with open(statsfile, 'r') as outstatsfile:
            next(outstatsfile)
            for fileline in outstatsfile:
                filelinedata = fileline.strip().split()
                statsdata["JobID"].append(filelinedata[0])
                statsdata["Elapsed"].append(filelinedata[1])
                statsdata["AllocCPUS"].append(filelinedata[2])
                statsdata["AveCPU"].append(filelinedata[3])
                statsdata["ReqMem"].append(filelinedata[4])
                statsdata["MaxVMSize"].append(filelinedata[5])
                statsdata["MaxRSS"].append(filelinedata[6])
                statsdata["MaxDiskRead"].append(filelinedata[7])
                statsdata["MaxDiskWrite"].append(filelinedata[8])
    except IOError:
        print(f"Could not read combined out stats data file: {statsfile}")
    finally:
        return statsdata


def determine_avgs_medians(statlist):
    list_avg = statistics.mean(statlist)
    list_median = statistics.median(statlist)
    list_stdev = statistics.stdev(statlist)
    return [list_avg, list_median, list_stdev]


def elapsedtime_to_seconds(elapsevals):
    elapsed_seconds = []
    for elapsevalue in elapsevals:
        elapsedata = elapsevalue.split(":")
        hts = int(elapsedata[0]) * 3600
        mts = int(elapsedata[1]) * 60
        elapsed_seconds.append(hts + mts + int(elapsedata[2]))
    return elapsed_seconds


def time_avg_median_stdev(elapsedtimevals):
    time_avg = statistics.mean(elapsedtimevals)
    time_median = statistics.median(elapsedtimevals)
    time_stdev = statistics.stdev(elapsedtimevals)
    return [time_avg, time_median, time_stdev]


def main():
    jds_params = get_params()
    statsdata = read_stats(jds_params["infile"])
    elapse_secs = elapsedtime_to_seconds(statsdata["Elapsed"])
    timestats = time_avg_median_stdev(elapse_secs)

    print("Overall time statistics for " +jds_params["jobname"]+ ":")
    print(f"AVG time: {round(timestats[0], 3)}s | {round(timestats[0]/60, 3)}m | {round(timestats[0]/3600, 3)}h")
    print(f"MED time: {round(timestats[1], 3)}s | {round(timestats[1]/60, 3)}m | {round(timestats[1]/3600, 3)}h")
    print(f"SD time: {round(timestats[2], 3)}s | {round(timestats[2]/60, 3)}m | {round(timestats[2]/3600, 3)}h")


if __name__ == "__main__":
    main()
