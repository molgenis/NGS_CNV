def filter_ccrs_calls(infileloc, outfileloc, filterneutrals, filtersize=0):
    """Filter CallCopyRatioSegments file and write filtered data to an output file.

    Parameters
    ----------
    infileloc : str
        Path to CallCopyRatioSegments file
    outfileloc : str
        Path to write output file to
    filterneutrals : bool
        Whether to filter out neutral calls
    filtersize : int
        Minimum CNV sizes to keep
    """
    try:
        outfile = open(outfileloc, 'w')
        with open(infileloc, 'r') as infile:
            for fileline in infile:
                if not fileline.startswith(("@", "CONTIG")):
                    write_line = True
                    filelinedata = fileline.strip().split("\t")

                    # Check whether the line is a to filter 
                    if filterneutrals:
                        if filelinedata[5] == '0':
                            write_line = False

                    # Check whether the call is smaller than the filter size
                    callsize = int(filelinedata[2]) - int(filelinedata[1])
                    if callsize <= filtersize:
                        write_line = False

                    # Check whether to write the line to output file
                    if write_line:
                        outfile.write(fileline)
                else:
                    outfile.write(fileline)
    except IOError:
        print(f"Could not filter {infileloc}")
