def remove_nas(infileloc, outfileloc):
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
