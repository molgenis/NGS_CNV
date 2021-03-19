# USAGE:
# THIS SCRIPT LINKS GENE NAMES BACK TO THE REGIONS IN THE HC_target.bed FILE.
# FIRST MAKE AN ALT_SLICED_BED FILE USING MY `alt_slice_bed_file.py` TO MAKE A SLICED BED FILE WITH GENE NAMES RETAINED
import sys

# Read the alt_sliced_merged.bed file.
alt_sliced_bed = {}
try:
    with open(sys.argv[1], 'r') as altsliced:
        for altline in altsliced:
            altlinedata = altline.strip().split("\t")
            region_str = f"{altlinedata[0]}:{altlinedata[1]}-{altlinedata[2]}"
            if region_str not in alt_sliced_bed:
                alt_sliced_bed[region_str] = altlinedata[3]
except IOError:
    print(f"Could not open {sys.argv[1]}")


# Check that alt sliced data has been read. If so loink gene names to regions in HC_target.bed
if len(alt_sliced_bed) > 0:
    try:
        with open(sys.argv[2]) as hcbedfile:
            for hcline in hcbedfile:
                hclinedata = hcline.strip().split("\t")
                hcregion = f"{hclinedata[0]}:{hclinedata[1]}-{hclinedata[2]}"
                if hcregion in alt_sliced_bed:
                    print(f"{hclinedata[0]}\t{hclinedata[1]}\t{hclinedata[2]}\t{alt_sliced_bed[hcregion]}")
    except IOError:
        print(f"Could not open HC_target bed file {sys.argv[2]}")
