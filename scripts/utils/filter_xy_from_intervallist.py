def filter_intervallist(intervalfile, outputfile):
    """Filter a provided interval list file to remove X and Y chromosome entries.

    Parameters
    ----------
    intervalfile : str
        Path to the interval list to filter
    outputfile : str
        Path to write filtered interval list file to
    """
    try:
        infile = open(intervalfile, "r")
        outfile = open(outputfile, "w")

        # Write the line in the input file to the output file if it is not an entry for chr X or Y
        for fileline in infile:
            if not fileline.startswith(XY_CHROMS):
                outfile.write(fileline)

        infile.close()
        outfile.close()
    except IOError:
        print("Could process interval list file. :(")
