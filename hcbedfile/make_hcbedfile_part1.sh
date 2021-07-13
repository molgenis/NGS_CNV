#Variables:
#$hcdir: Location of the main folder to use for creating the High Confident BED file
#$bedfile: Location of the BED file to start with when creating the HC BED file
#$slicedbedfile: Location and name of the sliced BED file, produced by slice_bed_file.py
#$umcuscriptsdir: Location of the UMCU Python scripts
#$umcuvenvdir: Location of the UMCU Python Virtual Env
#$generatesambamba: Directory containing the generate_sambamba.py script

echo "[...ACTIVATING UMCU PYTHON ENVIRONMENT...]"
source "${umcuvenvdir}"bin/activate

echo "[...SLICING BED FILE INTO INTERVALS OF 300...]"
python "${umcuscriptsdir}"slice_bed_file.py "${bedfile}" 300 > "${slicedbedfile}"

echo "[...MAKING POPULATION DIRECTORIES...]"
mkdir "${hcdir}"population1/male
mkdir "${hcdir}"population1/female
mkdir "${hcdir}"population2/male
mkdir "${hcdir}"population2/female

echo "[...PLACING MALE SAMPLES IN POPULATION1 DIRECTORY...]"
#<LINK M1>

echo "[...PLACING FEMALE SAMPLES IN POPULATION1 DIRECTORY...]"
#<LINK F1>

echo "[...PLACING MALE SAMPLES IN POPULATION2 DIRECTORY...]"
#<LINK M2>

echo "[...PLACING FEMALE SAMPLES IN POPULATION2 DIRECTORY...]"
#<LINK F2>

echo "[...GENERATING SAMBAMBA JOBS...]"
mkdir "${hcdir}"sambamba_jobs
python "${generatesambamba}"generate_sambamba.py "${hcdir}"population1/male "${hcdir}"sambamba_jobs m1
python "${generatesambamba}"generate_sambamba.py "${hcdir}"population1/female "${hcdir}"sambamba_jobs f1
python "${generatesambamba}"generate_sambamba.py "${hcdir}"population2/male "${hcdir}"sambamba_jobs m2
python "${generatesambamba}"generate_sambamba.py "${hcdir}"population2/female "${hcdir}"sambamba_jobs f2

echo "[...SUBMIT THE SAMBAMBA JOBS...]"
sh "${hcdir}"sambamba_jobs/submit_m1.sh
sh "${hcdir}"sambamba_jobs/submit_f1.sh
sh "${hcdir}"sambamba_jobs/submit_m2.sh
sh "${hcdir}"sambamba_jobs/submit_f2.sh

echo "[...DEACTIVATING UMCU PYTHON ENVIRONEMT...]"
deactivate
