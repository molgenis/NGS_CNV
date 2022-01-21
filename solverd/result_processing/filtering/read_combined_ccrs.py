from ccrscall import CcrsCall

def read_combined_ccrs(ccrsfileloc):
    """Read the combined CCRS data."""
    ccrsdata = {}
    ccrs_header_line = ""
    try:
        with open(ccrsfileloc, 'r') as ccrsfile:
            ccrs_header_line = next(ccrsfile)
            for fileline in ccrsfile:
                filelinedata = fileline.strip().split("\t")

                # Add the sample name and chromosome as dictionary entries
                if filelinedata[0] not in ccrsdata:
                    ccrsdata[filelinedata[0]] = {}
                if filelinedata[1] not in ccrsdata[filelinedata[0]]:
                    ccrsdata[filelinedata[0]][filelinedata[1]] = []
                ccrs_call = CcrsCall(filelinedata[0], filelinedata[1], int(filelinedata[2]), int(filelinedata[3]), int(filelinedata[4]), filelinedata[5], float(filelinedata[6]))

                # Check if there is frequency annotation
                if len(filelinedata) >= 9:
                    ccrs_call.ccrs_occurrence = filelinedata[7]
                    ccrs_call.ccrs_frequency = float(filelinedata[8])

                # Check if there is group frequency annotation
                if len(filelinedata) >= 12:
                    ccrs_call.ccrs_callgroup_name = filelinedata[9]
                    ccrs_call.ccrs_callgroup_occurrence = filelinedata[10]
                    ccrs_call.ccrs_callgroup_frequency = float(filelinedata[11])

                # Check if there is conrad annotation
                if len(filelinedata) >= 14:
                    ccrs_call.conrad_occurrence = filelinedata[12]
                    ccrs_call.conrad_frequency = float(filelinedata[13])

                # Check if there is gnomad annotation
                if len(filelinedata) >= 15:
                    ccrs_call.gnomad_frequency = float(filelinedata[14])

                # Add the CCRS call to the CCRS data
                ccrsdata[filelinedata[0]][filelinedata[1]].append(ccrs_call)
    except IOError:
        print("Could not read combined CCRS file :(")
    finally:
        return [ccrs_header_line, ccrsdata]
