def get_total_calls(infileloc, onlyna=False):
    """Simply determine and return the number of GATK4 CNV calls in the file.

    Parameters
    ----------
    infileloc : str
        Path to input file to read
    onlyna : bool
        If True only count lines with NA array CNV, if False count all lines
    """
    calls_total = 0
    try:
        with open(infileloc, 'r') as infile:
            next(infile)
            for fileline in infile:
                if onlyna:
                    filelinedata = fileline.strip().split("\t")
                    if filelinedata[4] == "NA":
                        calls_total += 1
                else:
                    calls_total += 1
        print(f"Total calls: {calls_total}")
    except IOError:
        print(f"Could not read ")
