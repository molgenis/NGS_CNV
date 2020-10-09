def remove_nas(infileloc, outfileloc):
    """Remove not overlapping with an array CNV from a classifications file.

    Parameters
    ----------
    infileloc : str
        Path to classifications file
    outfileloc: str
        Path to write filtered classifications output file to
    """
    try:
        outfile = open(outfileloc, 'w')
        with open(infileloc, 'r') as infile
            outfile.write(next(infile))
            for fileline in infile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[4] != "NA":
                    outfile.write(fileline)
        outfile.close()
    except IOError:
        print("Could not remove NAs :(")
