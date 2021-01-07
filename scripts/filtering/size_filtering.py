def determine_header_field(headerfields, columnname):
    """Identifies the index for the column name in the header.

    Parameters
    ----------
    headerfields : list of str
        List of header column names
    columnname : str
        Column name to fetch index for

    Returns
    -------
    indexfield : None or int
        Column index, None if column name not present in header
    """
    indexfield = None
    if columnname in headerfields:
        indexfield = headerfields.index(columnname)
    return indexfield


def filter_by_size(inputfile, outfileloc, callcolumn, minimumsize):
    """Filter a file with CNV calls to only retain calls satisfying a set minimum size.

    Parameters
    ----------
    inputfile : str
        Path to classification file to filter
    outfileloc : str
        Path to write filtered output file to
    callcolumn : str
        Column name containing the call 
    minimumsize : int
        Minimum size calls need to have to be retained

    Returns
    -------
    wrote_file : bool
        True if output file has been written, False if not
    """
    try:
        wrote_file = False
        outfile = open(outfileloc, 'w')

        with open(inputfile, 'r') as infile:
            headerline = next(infile)
            index_field = determine_header_field(headerline.strip.split("\t"), callcolumn)

            if index_field:
                outfile.write(headerline)
                for fileline in infile:
                    filelinedata = fileline.strip().split("\t")

                    if filelinedata[index_file] >= minimumsize:
                        outfile.write(fileline)
        outfile.close()
    except IOError:
        print("Could not filter input file.")
    finally:
        return wrote_file
