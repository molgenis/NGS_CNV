import sys
import os

indirfiles = os.listdir(sys.argv[1])
bamfiles = [f"{sys.argv[1]}/{bamfile}" for bamfile in indirfiles if bamfile.endswith(".bam")]
jobnum=0

for bamfile in bamfiles:
    jobnum += 1
    outfilepath = f"{sys.argv[2]}/sambamba_{jobnum}.sh"
    try:
        outfile = open(outfilepath, 'w')
        outfile.write("#!/bin/bash\n")
        outfile.write(f"#SBATCH --job-name=sambamba_{jobnum}\n")
        outfile.write(f"#SBATCH --output=sambamba_{jobnum}.out\n")
        outfile.write(f"#SBATCH --error=sambamba_{jobnum}.err\n")
        outfile.write("#SBATCH --time=04:00:00\n")
        outfile.write("#SBATCH --cpus-per-task=1\n")
        outfile.write("#SBATCH --mem=8gb\n")
        outfile.write("#SBATCH --nodes=1\n")
        outfile.write("#SBATCH --open-mode=append\n")
        outfile.write("#SBATCH --export=NONE\n")
        outfile.write("#SBATCH --get-user-env=L\n\n")
        outfile.write("module load sambamba/0.7.0\n")
        outfile.write("module list\n\n")
        outfile.write(f"/groups/umcg-gdio/tmp01/umcg-mbeukers/sambamba_0.6.5/sambamba_v0.6.5 depth region {bamfile} -L /groups/umcg-gdio/tmp01/umcg-mbeukers/highconfident_bedfile/sliced_captured.merged.bed -m -q 10 -F \"mapping_quality >= 20 and not duplicate and not failed_quality_control and not secondary_alignment\" -t 6 -o {bamfile}_coverage\n")
        outfile.close()
    except IOError:
        print(f"Could not write {outfile}")
