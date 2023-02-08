# NOTICE: We are updating this description, there might be some difference between the program and the text. Will fix this soon. (02/07/2023)

# LongAmpseq
The analysis for "CRISPR/Cas9 gene-editing of HSPCs from SCD patients has unintended consequences". 
The pipeline for illumina based large deletion profiling (LongAmp-seq) and Nano-pore based long-read analysis are stored in the following folders:  

-update 7/06/2022  

## Repository structure
**Folder: Illumina**  

  
**Folder: Nanopore**    
The scripts for processing Nanopore data (by Yilei Fu @ Treagen Lab).

## Environment set up
We recommend users to utilize virtual environment control tools such as conda to set up python environment.
The current version is on python 3, with the commands as follows for downloading required packages using conda:  
```
conda create -n longamp python=3
conda activate longamp
conda install -c bioconda bwa samtools seqtk bedtools bcftools
conda install numpy scikit-learn pandas scipy
conda install -c conda-forge matplotlib autoconf
```  

After setting up the environment, you can download the code for LongAmp-seq through git using the command below 
or download the zip file and unzip directly. 
```
git clone https://github.com/baolab-rice/LongAmpseq.git
cd LongAmpseq/Illumina
``` 

Then you will be ready for processing LongAmp-seq data.

## LongAmp-seq analysis

Before you start, please move the demultiplexed fastq files to the Illumina folder.
(Note that the standard output fastq of bcl2fastq should end with R1_001.fastq or R2_001.fastq)  
Note that you will need the reference genome ready for bwa alignment, 
and the path to the reference geneome in .fasta format will be required for LongAmp-seq analysis. 

+ Raw data processing  
 
For processing raw illumina fastqs (using common cell types that can utilize the reference genome directly)
```
bash longamp_bwa_PCR.sh [chromosome] [start index of long-range PCR] [index of cut site] [index of cut site +1] [end index of long-range PCR] [directory to the reference genome]
```
For example:
```
bash longamp_bwa_PCR.sh chr11 5245453 5248229 5248230 5250918 ~/Desktop/genomes/hg19/hg19.fa
```
And the output in csv format (file name ended with "filtered_2+.csv") should look like:

| Chromosome | Start    | End      | Read ID                                      | Score | Strand |
|------------|----------|----------|----------------------------------------------|-------|--------|
| chr11      | 59136655 | 59136956 | M04808:132:000000000-CTP75:1:2107:19955:2557 | 60    | -      |
| chr11      | 59138619 | 59138657 | M04808:132:000000000-CTP75:1:2107:19955:2557 | 60    | -      |

For reporter cell lines, please use 
```
bash longamp_bwa_cellline.sh GFPBFP 1 1004 1005 9423 GFPBFP_9423bp_NEW_ref_cutsite_1004bp.fa
```

With standard pipeline, three folders will be created:  
raw_data folder: containing original fastq files.  
processing folder: containing all the intermediate files (including the "filteredsorted.bam" used for delly calling).  
output folder: containing output for LongAmp-seq analysis listed as below:
```
filtered_2+.fastq: the collection of filtered full-length reads that will split in alignments  
filtered_2+.csv: the alignment pattern of split reads
filtered_1.fastq: the collection of filtered full-length reads that will NOT split in alignments

largedel_output.csv: all the read IDs containing large deletions
largedel_group.csv: large deletion patterns grouped by their start positions and lengths.bash longamp_bwa_cellline.sh GFPBFP 1 1004 1005 9423 GFPBFP_9423bp_NEW_ref_cutsite_1004bp.fa

```

+ Large deletion Profile generating  

The raw data processing section will call the correlated script and run it for you. 
For processing this step by yourselves, users can run the command as below:
```
# for general use
python bedfile.py [the filtered_2+.csv file] [output1: largedel_output.csv] [output2: largedel_group.csv] [length of amplicon]

# or for using reporter cell line
# Step 1: get the read number that spanning the cut site without large deletion:
echo $(cat [your file end with 30_filtered_1.fastq]|wc -l)/4|bc
# Step 2: run the following command
python bedfile_cellline.py [the filtered_2+.csv file] [output1: largedel_output.csv] [output2: largedel_group.csv] [length of amplicon, default=9423] [The number you got in Step 1]
```
After running bedfile.py commands you should be able to generate two csv files:  
largedel_output.csv: all the read IDs containing large deletions
largedel_group.csv: large deletion patterns grouped by their start positions and lengths.


```
# This command will process all the filtered and sorted bam files in the same directory.
bash longamp_delly_calling.sh [directory to the reference genome]
```
Therefore you will get the output with file name end with filteredsorted.csv,
which includes all the variants that called by delly.
This file will be analyzed for visualization in the next step.

+ Visualization  

The raw data processing section will call the correlated script and run it for you. 
For processing this step by yourselves, users can run the command as below:
```
python longampfigures_distribution.py [the largedel_output.csv file]
```
You can also specify the chromosomal index for the script:
```
python longampfigures_distribution.py [the largedel_output.csv file] [chromosome] [cut site] [length of long-range PCR]
```
For example:
```
python longampfigures_distribution.py [the largedel_output.csv file] chr11 5248229 5465
```

---------------------------------------

Please contact Yidan Pan (yidan.pan@rice.edu) and Mingming Cao (mc131@rice.edu) if you have any questions.
