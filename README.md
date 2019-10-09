# Macaques_SIV

# Steps to assess mutation frequencies in SIV-TB coinfected Macaques from R21 Study

#First after downloading all fastq files, move individual fastq files to the appropriate directory from individual subdirectories:

cp -r **/*.fastq.gz ./  or
mv -r **/*.fastq.gz ./  


# create the index for the reference genome before running the bashscripts

bwa index -p SIV -a bwtsw Data/SIV_Env.fasta



