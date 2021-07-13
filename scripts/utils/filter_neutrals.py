def filter_neutrals(infileloc, outfileloc):
    try:
        outfile = open(outfile, 'w')
        with open(infileloc, 'r') as infile:
        outfile.close()
    except:
        print(f"Could not filter neutrals from {infileloc}")
