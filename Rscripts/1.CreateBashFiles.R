#create the bash files to run bbmap and bwa
# read the template command text file:
cmmd<-readLines("Data/template/Bashxx.sh")

#choose the fastq files to be prrocessed
fq<-list.files("Data/AmbroseB670-2/", pattern="fastq.gz$") 
dir.create("Data/Bashscripts/")

#create vector of odd numbers:
n<-seq(1, by = 2, len = (length(fq)/2))
fq2<-fq[n]
for (i in 1:length(fq2)){
  #choose the paired reads fastq files
  fa1<-fq2[i]
  fa2<-gsub(pattern="R1",replace="R2",x=fa1)
  fname<-paste0("M",substr(fa1,start=1,stop=2))
  new<-gsub(pattern="10_S10_L001_R1_001.fastq", replace=paste0(fa1),x=cmmd)
  new<-gsub(pattern="10_S10_L001_R2_001.fastq", replace=paste0(fa2),x=new)
  new<-gsub(pattern="M10",replace=paste0(fname),x=new)
  writeLines(new, con=paste0("Data/Bashscripts/",fname,".sh"))

}

