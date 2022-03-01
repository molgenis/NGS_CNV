#!/usr/bin/env python
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


def filter_na(inputfile, outfileloc, colname):
    """Filter out calls without an overlapping Array CNV

    Parameters
    ----------
    inputfile : str
        CNV classification file to filter NA from
    outfileloc : str
        Path to write output file to
    colname : str
        Name of the column to use for filtering

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
            index_field = determine_header_field(headerline.strip().split("\t"), colname)

            if index_field:
                outfile.write(headerline)
                for fileline in infile:
                    filelinedata = fileline.strip().split("\t")
                    if filelinedata[index_field] != "NA":
                        outfile.write(fileline)
        outfile.close()
        wrote_file = True
    except IOError:
        print("Could not filter classifications file.")
    finally:
        return wrote_file
