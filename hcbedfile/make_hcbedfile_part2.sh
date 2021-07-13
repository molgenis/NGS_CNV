#Variables:
#$hcdir: Location of the main folder to use for creating the High Confident BED file
#$slicedbedfile: Location and name of the sliced BED file, produced by slice_bed_file.py
#$umcuscriptsdir: Location of the UMCU Python scripts
#$umcuvenv: Location of the UMCU Python Virtual Env

echo "[...ACTIVATING UMCU PYTHON ENVIRONMENT...]"
source "${umcuvenvdir}"bin/activate

echo "[...CALCULATE STATISTICS FOR EACH POPULATION...]"
python "${umcuscriptsdir}"filter_probe_file.py "${hcdir}"population1/male > "${hcdir}"population1_male_output_all
python "${umcuscriptsdir}"filter_probe_file.py "${hcdir}"population1/female > "${hcdir}"population1_female_output_all
python "${umcuscriptsdir}"filter_probe_file.py "${hcdir}"population2/male > "${hcdir}"population2_male_output_all
python "${umcuscriptsdir}"filter_probe_file.py "${hcdir}"population2/female > "${hcdir}"population2_female_output_all

echo "[...CALCULTE THE NUMBER OF AUTOSOMAL+X AND Y TARGETS IN SLICED BED FILE...]"
autoChrX=$(cat "${slicedbedfile}" | awk '($1 != "Y")' | wc -l)
chrY=$(cat "${slicedbedfile}" | awk '($1 == "Y")' | wc -l)

echo "[...CALCULATE THE NUMBER OF TARGETS TO REMAIN...]"
autoChrXremain=$(("${autoChrX}"*95/100))
chrYremain=$(("${chrY}"*33/100))

echo "[...FILTER MALE POPULATIONS...]"
cat "${hcdir}"population1_female_output_all | sed 's/inf/99999/g' | awk '($1 != "Y")' | sort -nk6 | head -n "${autoChrXremain}" | sed 's/X/999999999/g' | sort -nk1 -nk2 | sed 's/999999999/X/g' | awk '{OFS="\t"; print $1,$2,$3,$4"_"$5"_"$6}' > "${hcdir}"population1_female_auto_chrX
cat "${hcdir}"population2_female_output_all | sed 's/inf/99999/g' | awk '($1 != "Y")' | sort -nk6 | head -n "${autoChrXremain}" | sed 's/X/999999999/g' | sort -nk1 -nk2 | sed 's/999999999/X/g' | awk '{OFS="\t"; print $1,$2,$3,$4"_"$5"_"$6}' > "${hcdir}"population2_female_auto_chrX

echo "[...FILTER FEMALE POPULATIONS...]"
cat "${hcdir}"population1_male_output_all | sed 's/inf/99999/g' | awk '($1 == "Y")' | sort -nk6 | head -n "${chrYremain}" | sed 's/X/999999999/g' | sort -nk1 -nk2 | sed 's/999999999/X/g' | awk '{OFS="\t"; print $1,$2,$3,$4"_"$5"_"$6}' > "${hcdir}"population1_male_chrY
cat "${hcdir}"population2_male_output_all | sed 's/inf/99999/g' | awk '($1 == "Y")' | sort -nk6 | head -n "${chrYremain}" | sed 's/X/999999999/g' | sort -nk1 -nk2 | sed 's/999999999/X/g' | awk '{OFS="\t"; print $1,$2,$3,$4"_"$5"_"$6}' > "${hcdir}"population2_male_chrY

echo "[...MERGE MALE AND FEMALE POPULATION INTO ONE...]"
cat "${hcdir}"population1_female_auto_chrX "${hcdir}"population1_male_chrY > "${hcdir}"HC_population1
cat "${hcdir}"population2_female_auto_chrX "${hcdir}"population2_male_chrY > "${hcdir}"HC_population2

echo "[...CHECK THAT BOTH MERGED POPULATIONS HAVE THE SAME NUMBER OF ENTRIES...]"
firstPopEntries=$(wc -l "${hcdir}"HC_population1)
secondPopEntries=$(wc -l "${hcdir}"HC_population2)

echo "[...CALCULATE THE NUBMER OF OVERLAPPING ENTRIES BETWEEN THE TWO POPULATIONS...]"
popOverlapEntries=$(cat "${hcdir}"HC_populatie1 "${hcdir}"HC_populatie2 | cut -f1,2,3 | sort | uniq -c | awk '($1==2)' | wc -l)

echo "[...CALCULATE THE PERCENTAGE OVERLAP BETWEEN THE TWO POPULATIONS...]"
popOverlapThreshold=$(bc <<< 'scale=2; 99/100*100')
popOverlapLooseThreshold=$(bc <<< 'scale=2; 98/100*100')
popOverlapPercentage=$(bc <<< "scale=2; ${popOverlapEntries}/(${autoChrXremain}+${chrYremain})*100")
#popOverlapPercentage=$(("${popOverlapEntries}" / ("${autoChrXremain}" + "${chrYremain}") * 100))

echo "[...CHECKING WHETHER THE POPULATION OVERLAP SATISFIES THE STRICT OR LOOSE THRESHOLD...]"
if [ $(bc <<< "${popOverlapPercentage}>=${popOverlapThreshold}") -eq 1 ]
then
	echo "Population overlap satisfies 99% threshold so we can make the High Confident BED file"
	echo "[...CREATE THE HIGH CONFIDENT BED FILE...]"
	cat "${hcdir}"HC_populatie1 "${hcdir}"HC_populatie2 | cut -f1,2,3 | sort | uniq -c | awk '($1==2)' | sed 's/ \+/\t/g'  | cut -f 3,4,5 | sed 's/X/999999999/g'| sed 's/Y/9999999999/g' | sort -nk1 -nk2 | sed 's/9999999999/Y/g' | sed 's/999999999/X/g' > "${hcdir}"HC_target.bed
elif [ $(bc <<< "${popOverlapPercentage}>=${popOverlapLooseThreshold}") -eq 1 ]
then
	echo "Population overlap was lower than 99% but higher than 98%, therefore we will still make the High Confident BED file"
	echo "[...CREATE THE HIGH CONFIDENT BED FILE...]"
	cat "${hcdir}"HC_populatie1 "${hcdir}"HC_populatie2 | cut -f1,2,3 | sort | uniq -c | awk '($1==2)' | sed 's/ \+/\t/g'  | cut -f 3,4,5 | sed 's/X/999999999/g'| sed 's/Y/9999999999/g' | sort -nk1 -nk2 | sed 's/9999999999/Y/g' | sed 's/999999999/X/g' > "${hcdir}"HC_target.bed
else
	echo "Population overlap too small to properly make a High Confident BED file :("
fi

echo "[...DEACTIVATING UMCU PYTHON ENVIRONMENT...]"
deactivate
