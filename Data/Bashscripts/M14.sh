#!/bin/bash


mkdir Output/M14
#0 adapter trimming
bbduk.sh in1=Data/AmbroseB670-2/14_S14_L001_R1_001.fastq.gz in2=Data/AmbroseB670-2/14_S14_L001_R2_001.fastq.gz  out=Output/M14/M14_adp.trimmed.fastq ref=/Users/kahotisthammer/programs/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 stats=Output/M14/stats30_0.txt

#1 Trim reads at both ends at Score<30
bbduk.sh in=Output/M14/M14_adp.trimmed.fastq out=Output/M14/M14_trimmed.q30.fastq qtrim=rl trimq=30 stats=Output/M14/stats30_1.txt

#2. Kmer filtering
bbduk.sh in=Output/M14/M14_trimmed.q30.fastq out=Output/M14/M14_unmatched.q30.fq ref=/Users/kahotisthammer/programs/bbmap/resources/phix174_ill.ref.fa k=31 hdist=1 stats=Output/M14/stats30_2.txt

#3.Remove reads with ave score <30 (Phred score =20 99% accuracy 1% chance of mistake)
bbduk.sh in=Output/M14/M14_unmatched.q30.fq out=Output/M14/M14_clean.q30.fq maq=30 stats=Output/M14/stats30_3.txt

#4. Align the file using bwa to the reference 

bwa mem -t 4 -k 15 -a SIV Output/M14/M14_clean.q30.fq  > Output/M14/M14_BWAmapped.sam

#5. convert sam to bam
samtools view -S -b Output/M14/M14_BWAmapped.sam > Output/M14/M14_BWAmapped.bam


#6. sort the bam file
samtools sort  Output/M14/M14_BWAmapped.bam -o  Output/bam/M14_BWA_sort.bam

#7. index the bam file
samtools index  Output/bam/M14_BWA_sort.bam  Output/bam/M14_BWA.sort_bam.bai

rm Output/M14/M14_adp.trimmed.fastq Output/M14/M14_trimmed.q30.fastq Output/M14/M14_unmatched.q30.fq Output/M14/M14_BWAmapped.sam

