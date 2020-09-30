def parse_regions(genomicregions, paddingtoadd):
    """Parse one or more genomic regions and return the parsed data.

    Parameters
    ----------
    genomicregions : str or list of str
        Genomic region(s) to parse.

    Returns
    -------
    region_data : dict
        Parsed genomic region data.
    """
    region_data = {}

    # Check if one or multiple genomic regions have been specified.
    if type(genomicregions) is not list:
        genomicregions = [genomicregions]

    # Iterate over the regions and save the parsed data.
    for genregion in genomicregions:
        genregion = parse_region(genregion, paddingtoadd)
        
        # Add the region as a range to the region data dict.
        if genregion[0] not in region_data:
            region_data[genregion[0]] = []
        region_data[genregion[0]].append(range(genregion[1], genregion[2]))
    return region_data


def parse_region(genomeregion, paddingtoadd):
    """Parse a genomic region and return the info as separate fields.

    Genomic regions are expected to be in the format chr:start-end.
    The individual parts returned are chrom, start and end.
    Padding is added to the start (via subtraction) and end (via addition) positions.

    Parameters
    ----------
    genomeregion : str
        Genomic region in the format chr:start-end
    paddingtoadd : int
        Padding to add to start and end
    """
    region_split = []
    chrom_pos = genomeregion.split(":")
    region_split.append(chrom_pos[0])    # Add the genome name
    region_split.append(int(chrom_pos[1].split("-")[0]) - paddingtoadd)    # Add the genomic region start - padding
    region_split.append(int(chrom_pos[1].split("-")[1]) + paddingtoadd)    # Add the genomic region end + padding
    return region_split


def extract_intervals(inputfile, genomicregions):
    """Extract and return interval data overlapping with a genomic region.

    Parameters
    ----------
    inputfile : str
        Path to TSV/SEG file containing intervals to extract.
    genomicregions : 
        Genomic regions to use for selecting intervals.

    Returns
    -------
    extracted_intervals : dict
        Intervals from the input file that overlap with the provided genomic regions.
    """
    extracted_intervals = {}
    try:
        with open(inputfile, "r") as infile:
            for fileline in infile:
                if not fileline.startswith(("@", "CONTIG")):
                    filelinedata = fileline.strip().split("\t")
                    overlap_region = get_interval_overlap(filelinedata[0], int(filelinedata[1]), int(filelinedata[2]), genomicregions)
                    
                    # Check if the interval overlaps with a genomic region. If so keep it.
                    if overlap_region:
                        if overlap_region not in extracted_intervals:
                            extracted_intervals[overlap_region] = []
                        extracted_intervals[overlap_region].append(fileline.strip())
    except IOError:
        print(f"Could not read input file: {inputfile}")
    finally:
        return extracted_intervals


def get_interval_overlap(intervalchrom, intervalstart, intervalend, genomicregions):
    """Return the genomic region the interval overlaps with. Return None if there is no overlap.

    Parameters
    ----------
    intervalchrom : str
        Interval chromosome name.
    intervalstart : int
        Interval start position.
    intervalend : int
        Interval end position.
    genomicregions : dict
        Genomic regions to check for overlap.

    Returns
    -------
    str or None
        Overlapping genomic region as chr:start-end, None if there is no overlap.
    """
    if intervalchrom in genomicregions:
        genomicranges = genomicregions[intervalchrom]
        for genomicrange in genomicranges:
            if intervalstart in genomicrange or intervalend in genomicrange:
                return f"{intervalchrom}:{genomicrange[0]}-{genomicrange[-1]}"
    return None


def make_interval_label(labelchrom, labelstart, labelend):
    """Create and return the genomic interval data as a string label.
    
    Parameters
    ----------
    labelchrom : str
        Chromosome name the interval is located on.
    labelstart : str
        Start position of the interval.
    labelend : str
        End position of the interval.
    """
    return f"{labelchrom}:{labelstart}-{labelend}"
