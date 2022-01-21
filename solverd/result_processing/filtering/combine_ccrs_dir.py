import sys
import os
from ccrscall import CcrsCall


def read_combined_ccrs(ccrsfileloc):
    """Read the combined CCRS data."""
    ccrsdata = {}
    try:
        with open(ccrsfileloc, 'r') as ccrsfile:
            global_ccrs_header_line = next(ccrsfile)
            for fileline in ccrsfile:
                filelinedata = fileline.strip().split("\t")

                # Add the sample name and chromosome as dictionary entries
                if filelinedata[0] not in ccrsdata:
                    ccrsdata[filelinedata[0]] = {}
                if filelinedata[1] not in ccrsdata[filelinedata[0]]:
                    ccrsdata[filelinedata[0]][filelinedata[1]] = []
                ccrs_call = CcrsCall(filelinedata[0], filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), int(filelinedata[4]), filelinedata[5], float(filelinedata[6]))

                # Check if the there is frequency annotation
                if len(filelinedata) >= 9:
                    ccrs_call.ccrs_occurrence = filelinedata[7]
                    ccrs_call.ccrs_frequency = float(filelinedata[8])

                # Check if there is conrad annotation
                if len(filelinedata) >= 11:
                    ccrs_call.conrad_occurrence = filelinedata[9]
                    ccrs_call.conrad_frequency = float(filelinedata[10])

                # Check if there is gnomad annotation
                if len(filelinedata) >= 12:
                    ccrs_call.gnomad_frequency = float(filelinedata[11])

                # Add the CCRS call to the CCRS data
                ccrsdata[filelinedata[0]][filelinedata[1]].append(ccrs_call)
    except IOError:
        print("Could not read combined CCRS file :(")
    finally:
        return ccrsdata


def write_combined_ccrs(outfileloc, ccrsdata):
    file_written = False
    sample_names = list(ccrsdata.keys())
    sample_names.sort()
    try:
        with open(outfileloc, 'w') as outfile:
            for samplename in sample_names:
                for chromname in ccrsdata[samplename]:
                    for ccrscall in ccrsdata[samplename][chromname]:
                        outfile.write(ccrscall.to_ccrs_file_line())
        file_written = True
    except IOError:
        print("Could not write combined CCRS file")
    finally:
        return file_written


ccrs_files = [x for x in os.listdir(sys.argv[1]) if x.endswith(".seg")]
all_ccrs_data = {}
for ccrsfile in ccrs_files:
    all_ccrs_data.update(read_combined_ccrs(f"{sys.argv[1]}{ccrsfile}"))
aap = write_combined_ccrs(sys.argv[2], all_ccrs_data)
print(aap)



# sys.argv[1]: directory with ccrs files
# sys.argv[2]: path to output file
