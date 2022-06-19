#!/bin/bash

start=${1:-5245453}
end=${2:-5250918}
filename=${3}

grep -v '^@' ${filename} | awk -v num1=${start} -v num2=${end} '{if ($1 ~ /^@/ && $3!="chr11" || $4<num1 || $4>num2) print}' > ${filename/.sam/_LI.sam}
echo -e "ReadID\tChr\tStart" > ${filename/.sam/_hit.txt}
awk -v OFS='\t' '{if ($3!="*") print $1"\t"$3"\t"$4}' ${filename/.sam/_LI.sam} | sort -n -k3 | sort -V -k2 | awk '!seen[$3]++' >> ${filename/.sam/_hit.txt} 
