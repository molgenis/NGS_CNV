#!/usr/bin/env python
def collect_allelic_log2_values(allelic_data, interval_data, outfilepath):
    """Collect log2 values for the allelic positions and write results to an output file.
    
    Parameters
    ----------
    allelic_data : dict
        Read allelic data.
    interval_data : dict
        Read interval data.
    outfilepath : str
        Path to write output file to.
    """
    try:
        with open(outfilepath, "w") as outfile:
            outfile.write("CONTIG\tPOSITION\tLOG2_COPY_RATIO\tOVERLAP\tLABEL\n")
            
            # Start collecting the log2 copy ratio value for the allelic positions.
            for allelic_line in allelic_data:
                allelic_linedata = allelic_data[allelic_line]
                log2value = get_log2_value(allelic_linedata[0], int(allelic_linedata[1]), interval_data)
                
                if log2value:
                    outfile.write(f"{allelic_linedata[0]}\t{allelic_linedata[1]}\t{log2value[3]}\t{log2value[4]}\t{log2value[5]}\n")
                else:
                    outfile.write(f"{allelic_linedata[0]}\t{allelic_linedata[1]}\tNA\tNA\tNA\n")
    except IOError:
        print("Could not open ")


def get_log2_value(allelic_chrom, allelic_pos, interval_data):
    """Return the log2 copy ratio value for the allelic position.
    
    Parameters
    ----------
    allelic_chrom : str
        Chromosome name of the SNP allele.
    allelic_pos : int
        Position of the SNP allele.
    interval_data : dict
        Read intervals.
    
    Returns
    -------
    list
        Data of the overlapping interval.
    """
    # Start iterating over the intervals.
    for intervalkey in interval_data:
        intervalkey_data = intervalkey.split("_")
        intervalkey_chrom = intervalkey_data[0]
        intervalkey_pos = int(intervalkey_data[1])
        
        # Check whether the allelic position overlaps with the interval.
        if allelic_chrom == intervalkey_chrom:
            if allelic_pos == intervalkey_pos:
                return interval_data[intervalkey]
            elif allelic_pos > intervalkey_pos:
                if allelic_pos <= int(interval_data[intervalkey][2]):
                    return interval_data[intervalkey]
