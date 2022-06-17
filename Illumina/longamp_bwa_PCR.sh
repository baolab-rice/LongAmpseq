#!/bin/bash

chr_num=${1:-chr11}
pcr_start=${2:-5245453}
cut_site_left=${3:-5248229}
cut_site_right=${4:-5248230}
pcr_end=${5:-5250918}
ref_genome=${6:-/home/yp11/Desktop/genomes/hg19/hg19.fa}
##########################
# put the script in the same directory as demultiplexed fastq files.

# end with fastq
# start with merging
# use full hg19
# use q30 only for the initial alignment and filtering (same for 0915)
# removed the duplication deduct part

# Yidan Pan @Baolab
# Ming Cao @Baolab
##########################
  # update 06/16/2022

echo 'LongAmp-seq analysis' > log.txt

for file in *R1_001.fastq
do
  flash -M600 ${file} ${file/R1/R2} 2>&1 | tee flash.log # read length = 300 so we use 300*2 as maximum merged read length 
  mv ./out.extendedFrags.fastq ${file/.fastq/_merged.fastq}
done

echo "The total pairs found by FLASH" >> log.txt
grep "Total pairs:" flash.log | awk '{print $4}' >> log.txt
echo "The combined pairs done by FLASH" >> log.txt
grep "Combined pairs:" flash.log | awk '{print $4}' >> log.txt

for file in *merged.fastq
do
# init
  mkdir ${file/merged.fastq/longamp}
  mkdir ${file/merged.fastq/longamp}/raw_data
  mkdir ${file/merged.fastq/longamp}/processing
  mkdir ${file/merged.fastq/longamp}/output

#filtering
  bwa mem -A2 -E1 ${ref_genome} ${file} > ${file/.fastq/.sam}
  samtools view -S -b -q 30 ${file/.fastq/.sam} | samtools sort -o ${file/.fastq/_sorted.bam}
  samtools index ${file/.fastq/_sorted.bam}
  samtools view -b ${file/.fastq/_sorted.bam} ${chr_num}:${pcr_start}-${cut_site_left} > ${file/.fastq/_sortedleft.bam}
  samtools index ${file/.fastq/_sortedleft.bam}
  samtools view -b ${file/.fastq/_sorted.bam} ${chr_num}:${cut_site_right}-${pcr_end} > ${file/.fastq/_sortedright.bam}
  samtools index ${file/.fastq/_sortedright.bam}
  samtools view -F 4 ${file/.fastq/_sortedleft.bam} | cut -f1 | sort -u > ${file/.fastq/_leftID.txt}
  samtools view -F 4 ${file/.fastq/_sortedright.bam} | cut -f1 | sort -u > ${file/.fastq/_rightID.txt}
  comm -12 ${file/.fastq/_leftID.txt} ${file/.fastq/_rightID.txt} > ${file/.fastq/_bothID.txt}
  seqtk subseq ${file} ${file/.fastq/_bothID.txt} > ${file/.fastq/30_filtered.fastq}
#now we have filtered fastq

  bwa mem -A2 -E1 ${ref_genome} ${file/.fastq/30_filtered.fastq} >${file/.fastq/30_filtered.sam}
  samtools view -S -b ${file/.fastq/30_filtered.sam} -o ${file/.fastq/30_filtered.bam}
  samtools sort ${file/.fastq/30_filtered.bam} -o ${file/.fastq/30_filteredsorted.bam}
  bedtools bamtobed -i ${file/.fastq/30_filteredsorted.bam} > ${file/.fastq/30_filtered.bed}
#get the filtered, non-deducted reads

#extract the reads that have supplementary alignments
  samtools view -F 4 ${file/.fastq/30_filtered.bam} | cut -f1 | uniq -d > ${file/.fastq/filtered_2+alignID.txt}
  samtools view -F 4 ${file/.fastq/30_filtered.bam} | cut -f1 | uniq -u > ${file/.fastq/filtered_1alignID.txt}
  seqtk subseq ${file/.fastq/30_filtered.fastq} ${file/.fastq/filtered_2+alignID.txt} > ${file/.fastq/30_filtered_2+.fastq}
  seqtk seq -a ${file/.fastq/30_filtered_2+.fastq} > ${file/.fastq/30_filtered_2+.fasta}
  seqtk subseq ${file/.fastq/30_filtered.fastq} ${file/.fastq/filtered_1alignID.txt} > ${file/.fastq/30_filtered_1.fastq}

  bwa mem -A2 -E1 ${ref_genome} ${file/.fastq/30_filtered_2+.fastq} >${file/.fastq/filtered_2+.sam}
  samtools view -S -b ${file/.fastq/filtered_2+.sam} -o ${file/.fastq/filtered_2+.bam}
  bedtools bamtobed -i ${file/.fastq/filtered_2+.bam} > ${file/.fastq/filtered_2+.bed}

  #convert to csv
  cat ${file/.fastq/filtered_2+.bed} | tr "\\t" "," > ${file/.fastq/filtered_2+.csv}

  amp_length=$((pcr_end-pcr_start))
  python bedfile.py ${file/.fastq/filtered_2+.csv} ${file/.fastq/largedel_output.csv} ${file/.fastq/largedel_group.csv} ${amp_length}
  python add_seq.py ${file/.fastq/largedel_output.csv} ${file/.fastq/30_filtered_2+.fasta} 

  # clustering
  python3 clustering_LAS.py ${file/.fastq/largedel_output.csv} > ${file/.fastq/largedel_output_cluster.csv}
  ## filter cluster size (change to cluster size <0.1% read)
  #cluster_size=$(awk '{SUM+=$1} END{print SUM}' ${file/.fastq/largedel_output_cluster.csv} | xargs -I {} echo {}*0.001 | bc)
  ## filter with aligned read number
  read_num=$(echo $(cat ${file/.fastq/30_filtered.fastq}|wc -l)/4|bc)
  read_size=$(echo $read_num*0.0001 | bc)
  touch ${file/.fastq/largedel_output_cluster_t.csv}
  awk -v cluster=$read_size '{if ($1>cluster) print}' ${file/.fastq/largedel_output_cluster.csv} > ${file/.fastq/largedel_output_cluster_t.csv}
  cp ${file/.fastq/largedel_output_cluster_t.csv} ${file/.fastq/largedel_output_cluster.csv} 

  # visualization
  python3 distribution_LAS.py ${file/.fastq/largedel_output_cluster.csv} ${cut_site_left}
  sleep 2
  # append read numbers into the log file.
  # total aligned events:
  echo 'Total aligned read number after filtering:' >> log.txt
  echo ${file/_merged.fastq/} >> log.txt
  echo $(cat ${file/.fastq/30_filtered.fastq}|wc -l)/4|bc >> log.txt
  echo 'Read with unique alignment(use for crispresso2):' >>log.txt
  echo $(cat ${file/.fastq/30_filtered_1.fastq}|wc -l)/4|bc >> log.txt
  echo 'Read with more than one aligned segments:' >>log.txt
  echo $(cat ${file/.fastq/30_filtered_2+.fastq}|wc -l)/4|bc >> log.txt
  echo 'Read with large deletions:' >>log.txt
  echo $(cat ${file/.fastq/largedel_output.csv}|wc -l)-1|bc >> log.txt

  # file arrangement
  mv ${file/R1_001_merged.fastq/}R*_001.fastq ${file/merged.fastq/longamp}/raw_data/
  mv ${file} ${file/merged.fastq/longamp}/raw_data/

  mv ${file/.fastq/30_filtered_2+.fastq} ${file/merged.fastq/longamp}/output
  mv ${file/.fastq/30_filtered_2+.fasta} ${file/merged.fastq/longamp}/output
  mv ${file/.fastq/30_filtered_1.fastq} ${file/merged.fastq/longamp}/output
  mv ${file/.fastq/filtered_2+.csv} ${file/merged.fastq/longamp}/output
  mv ${file/.fastq/largedel_output.csv} ${file/merged.fastq/longamp}/output
  mv ${file/.fastq/largedel_group.csv} ${file/merged.fastq/longamp}/output
  #mv ${file/.fastq/30_filteredsorted.csv} ${file/merged.fastq/longamp}/output/${file/merged.fastq/_delly.csv}
  #mv ${file/.fastq/30_filteredsorted.svg} ${file/merged.fastq/longamp}/output/${file/merged.fastq/_delly.svg}
  mv ${file/.fastq/largedel_output_cluster.csv} ${file/merged.fastq/longamp}/output/${file/.fastq/largedel_output_cluster.csv}
  mv ${file/.fastq/largedel_output_cluster.svg} ${file/merged.fastq/longamp}/output/${file/.fastq/largedel_output_cluster.svg}    
  
  mv *${file/merged.fastq/}* ${file/merged.fastq/longamp}/processing 2>/dev/null
  mv log.txt ${file/merged.fastq/longamp}/output/${file/.fastq/_log.txt}
  mv flash.log ${file/merged.fastq/longamp}/output



done
