library(ggplot2)
library(reshape)
library(ggpubr)
library(ggthemes)
library(plotrix)
library(grid)

source("Rscripts/baseRscript2.R")

# read the files saved in Overview (unfiltered):
SIVFiles_overview<-list.files("Output/Overview/",pattern="overview.csv")

Overview_summary<-list()
for (i in 1:length(SIVFiles_overview)){ 
  overviews<-read.csv(paste0("Output/Overview/",SIVFiles_overview[i]),stringsAsFactors=FALSE, row.names = 1)
  Overview_summary[[i]]<-overviews
  names(Overview_summary)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=3)
}

# PrimerID prepped files
SIVFiles_overview499<-list.files("Output/Overview_499/",pattern="overview.csv")

Overview_summary499<-list()
for (i in 1:length(SIVFiles_overview499)){ 
  overviews<-read.csv(paste0("Output/Overview_499/",SIVFiles_overview499[i]),stringsAsFactors=FALSE, row.names = 1)
  Overview_summary499[[i]]<-overviews
  names(Overview_summary499)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=3)
}

pdf(file=paste0("Output/Comparison/TransitionMF_comparison.pdf"), width=15,height=18)
par(mfrow = c(4,2))
for (i in 1:length(Overview_summary)){
      S<-Overview_summary[[i]]
      S<-S[which(grepl(269, S$pos)):which(grepl(767, S$pos)),]
      S_499<-Overview_summary499[[i]]
      filename<-substr(names(Overview_summary)[i],start=1, stop=3)
      
      mf1<-S[,c("pos","freq.Ts")]
      mf2<-S_499[,c("pos","freq.Ts")]
      colnames(mf2)[2]<-"freq.Ts2"
      S.Ts<-merge(mf1,mf2,by="pos")
      write.csv(S.Ts, paste0("Output/Comparison/TsMutFreq_",filename,".csv"))
      
      #pdf(file=paste0("Output/Comparison/TransitionMF_",filename,".pdf"), width=7,height=4)
      plot(S.Ts$freq.Ts, pch=16,cex=0.8, col="#0077BB",xlab="Nucleotide position",ylab = "Transition mutation frequency", main= paste0("Mutation frequency comparison: ",filename))
      points(S.Ts$freq.Ts2, pch=1, col="#EE6677", cex=0.6)
      legend("topleft",legend=c("All", "PrimerID"), col=c("#0077BB","#EE6677"),pch=c(16,1),cex=0.8)
      
      mean1<-format(round(mean(S.Ts$freq.Ts,na.rm=T),7))
      mean2<-format(round(mean(S.Ts$freq.Ts2,na.rm=T),7))
      mtext(paste0("Ave.= ",mean1),side=1, at=450, line=2, col= "#0077BB", cex=0.8)
      mtext(paste0("Ave.= ",mean2),side=1, at=450, line=2.9, col= "#EE6677", cex=0.8)
                  
      
      
}
dev.off()


pdf(file="Output/Comparison/TsMF_difference.pdf", width=15,height=18)
par(mfrow = c(4,2))
for (i in 1:length(Overview_summary)){
      S<-Overview_summary[[i]]
      S<-S[which(grepl(269, S$pos)):which(grepl(767, S$pos)),]
      S_499<-Overview_summary499[[i]]
      filename<-substr(names(Overview_summary)[i],start=1, stop=3)
      mf1<-S[,c("pos","freq.Ts")]
      mf2<-S_499[,c("pos","freq.Ts")]
      colnames(mf2)[2]<-"freq.Ts2"
      S.Ts<-merge(mf1,mf2,by="pos")
      plot(S.Ts$freq.Ts-S.Ts$freq.Ts2, pch=16,cex=0.8, col="#0077BB",xlab="Nucleotide position", ylab = "Differences in mutation frequencies" ,main=paste0("Mutation frequency comparison: ",filename))
      ave<-mean(abs(S.Ts$freq.Ts-S.Ts$freq.Ts2), na.rm=T)
      ave<-format(round(ave,6))
      mtext(paste0("ave.diff=",ave), side=1, at =440, line=2.5, cex=0.9, col="#004488")
      }
dev.off()
      
      
