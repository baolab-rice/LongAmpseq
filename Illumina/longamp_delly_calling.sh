#!/bin/bash

ref_genome=${1:-/home/yp11/Desktop/genomes/hg19/hg19.fa}

for file in *filteredsorted.bam
  do samtools index ${file}
  chmod a+x delly
  ./delly call -o ${file/.bam/.bcf} -g ${ref_genome} ${file}

done

for file in *.bcf; do bcftools query -f'%CHROM\t%POS0\t%END\t%ID\n' ${file} > ${file/.bcf/.bed}; cat ${file/.bcf/.bed} | tr "\\t" "," > ${file/.bcf/.csv}; done

# for file in *.bcf; do cat ${file/.bcf/.bed} | tr "\\t" "," > ${file/.bcf/.csv}; done