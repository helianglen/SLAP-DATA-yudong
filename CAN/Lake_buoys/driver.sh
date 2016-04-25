
# 4/25/2016
# Buoy data were copy&pasted from Excel files: SLAP_Houser_Buoy_1.xlsx and _2.xlsx. 

for bid in buoy1 buoy2; do 

awk -F'\t' '{print $2, $3}' ${bid}.csv > ${bid}_temp.txt 
awk -F'\t' '{print $1}' ${bid}.csv > ${bid}_date.txt 
./compute_buoy_Tb ${bid}_temp.txt > ${bid}_Tb.txt 
pr -tm -s', ' ${bid}_date.txt ${bid}_Tb.txt > ${bid}_all_data.csv 

done






