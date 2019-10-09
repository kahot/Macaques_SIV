
SIVFiles_overview<-list.files("Output/Overview/",pattern="overview.csv")

Overview_summary<-list()
for (i in 1:length(SIVFiles_overview)){ 
  overview<-read.csv(paste0("Output/Overview/",SIVFiles_overview[i]),stringsAsFactors=FALSE, row.names = 1)
  Overview_summary[[i]]<-overview
  names(Overview_summary)[i]<-substr(paste(SIVFiles_overview[i]),start=1,stop=3)
}
#################
##############
# Filter out mut.freq over 0.2

FilteredOverview<-list()
for (i in 1:length(Overview_summary)){
  dat<-Overview_summary[[i]]
  filename<-names(Overview_summary)[i]
  
  #filter the bases at non-equilibrium (mut. freq over 0.2)
  dat$freq.Ts[dat$freq.Ts>=0.2] <-NA
  dat$freq.Ts.ref[dat$freq.Ts.ref>=0.2] <-NA
  dat$freq.transv1.ref[dat$freq.transv1.ref>=0.2] <-NA
  dat$freq.transv2.ref[dat$freq.transv2.ref>=0.2] <-NA
  dat$freq.transv.ref[dat$freq.transv.ref>=0.2] <-NA
  dat$freq.mutations.ref[dat$freq.mutations.ref>=0.2] <-NA
  
  FilteredOverview[[i]]<-dat
  write.csv(dat,paste0("Output/overview/",filename,"_filtered.overview2.csv"))
  names(FilteredOverview)[i]<-filename
}

#make the plots again 
source("Rscripts/MutationFreqSumFil.R")

#the results are saved in Output/MutFreq_fil/ directory

#1) Transistion mutations
TransFiles<-list.files("Output/MutFreq_fil/",pattern="Transition.csv")
TransFiles<-TransFiles[-6]

NofSamples<-length(FilteredOverview)

for (i in 1:length(TransFiles)){
  mdata<-read.csv(paste0("Output/MutFreq_fil/",TransFiles[i]),stringsAsFactors=FALSE,row.names=1)
  mdata<-mdata[1:NofSamples,]
  colnames(mdata)<-c("A","T","C","G")
  Dat<-melt(mdata)
  Dat<-Dat[!(is.na(Dat$value)),]
  
  filename<-sub(".csv$","",paste(TransFiles[i]))
  
  MFplot<-ggplot(Dat,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
    scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC")) + theme_classic()+
    guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
    theme(axis.text.y = element_text(size =10))+
    ggtitle(paste("Ave. freq of ",filename," mutation")) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
    theme(plot.title = element_text(hjust = 0.5))
  if (i==1|i==2){
    ggsave(filename=paste0("Output/MutFreq_fil/",filename,".pdf"),width=5, height=7, units='in',device='pdf', plot=MFplot)}
  else{
    ggsave(filename=paste0("Output/MutFreq_fil/",filename,".pdf"), width=9, height=7, units='in',device='pdf',plot=MFplot)
  }
  
}

########### Same Transition Mutation Figures but two plots side by side #########
#plots Syn vs NonSyn each other:
for (i in 5:6){
  mdata<-read.csv(paste0("Output/MutFreq_fil/",TransFiles[i]),stringsAsFactors=FALSE,row.names=1)
  mdata<-mdata[1:NofSamples,]
  colnames(mdata)<-c("A","T","C","G")
  Dat<-melt(mdata)
  Dat<-Dat[!(is.na(Dat$value)),]
  
  filename<-sub(".csv$","",paste(TransFiles[i]))
  
  MFplot<-ggplot(Dat,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
    scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC")) + theme_classic()+
    guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
    theme(axis.text.y = element_text(size =10))+
    ggtitle(paste(filename)) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
    theme(plot.title = element_text(hjust = 0.5))
  plotname<-paste0("plot_",i)
  
  assign(plotname,MFplot)
}
require(gridExtra)
pdf("Output/MutFreq_fil/SynvsNonsynTs.pdf",width=10, height=7)
grid.arrange(plot_5, plot_6, ncol=2)
dev.off()

# Plot Syn CpgMaking vs. NonCpG
for (i in c(2,4)){
  mdata<-read.csv(paste0("Output/MutFreq_fil/",TransFiles[i]),stringsAsFactors=FALSE,row.names=1)
  mdata<-mdata[1:NofSamples,]
  colnames(mdata)<-c("A","T","C","G")
  Dat<-melt(mdata)
  Dat<-Dat[!(is.na(Dat$value)),]
  
  filename<-sub(".csv$","",paste(TransFiles[i]))
  
  MFplot<-ggplot(Dat,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
    scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC")) + theme_classic()+
    guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
    theme(axis.text.y = element_text(size =10))+
    ggtitle(paste(filename)) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
    theme(plot.title = element_text(hjust = 0.5))
  plotname<-paste0("plot_",i)
  
  assign(plotname,MFplot)
}
pdf("Output/MutFreq_fil/CpG.vs.NonCpG.SynTs.pdf",width=10, height=7)
grid.arrange(plot_2, plot_4, ncol=2)
dev.off()


#################
#Plot summary of 2) Transversion mutations         
Tv1Files<-list.files("Output/MutFreq_fil/",pattern="Transversion_1.csv")
Tv2Files<-list.files("Output/MutFreq_fil/",pattern="Transversion_2.csv")

#CpGmaking Transversion plots
for (i in 1:2){
  mdata1<-read.csv(paste0("Output/MutFreq_fil/",Tv1Files[i]),stringsAsFactors=FALSE,row.names=1)
  mdata2<-read.csv(paste0("Output/MutFreq_fil/",Tv2Files[i]),stringsAsFactors=FALSE,row.names=1)
  mdata1<-mdata1[1:NofSamples,]
  mdata2<-mdata2[1:NofSamples,]
  mdata<-cbind(mdata1[,c(1,4)],mdata2[,c(2,3)])
  mdata<-mdata[,c(1,3,4,2)]
  colnames(mdata)<-c("A","T","C","G")
  dat<-melt(mdata)
  
  filename<-sub("_1.csv$","",paste(Tv1Files[i]))
  
  MFplot<-ggplot(dat,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
    scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC")) + theme_classic()+
    guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
    theme(axis.text.y = element_text(size =10))+
    ggtitle(paste0("Ave. freq of ",filename," mutation")) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename=paste0("Output/MutFreq_fil/",filename,".pdf"),width=9, height=7, units='in',device='pdf',plot=MFplot)
  
  dev.off()
}

#non-CpGmaking Transversion plots
for (i in c(3,4,5,7)){
  mdata1<-read.csv(paste0("Output/MutFreq_fil/",Tv1Files[i]),stringsAsFactors=FALSE,row.names=1)
  mdata2<-read.csv(paste0("Output/MutFreq_fil/",Tv2Files[i]),stringsAsFactors=FALSE,row.names=1)
  mdata1<-mdata1[1:NofSamples,]
  mdata2<-mdata2[1:NofSamples,]
  mdata<-mdata1+mdata2
  colnames(mdata)<-c("A","T","C","G")
  dat<-melt(mdata)
  
  filename<-sub("_1.csv$","",paste(Tv1Files[i]))
  
  MFplot<-ggplot(dat,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Nucleotide",y="Mutation frequency")+
    scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC")) + theme_classic()+
    guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
    theme(axis.text.y = element_text(size =10))+
    ggtitle(paste0("Ave. freq of ",filename," mutation")) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(filename=paste0("Output/MutFreq_fil/",filename,".pdf"),width=9, height=7, units='in',plot=MFplot)
  
  dev.off()
}        



#####
#Average all bases and plot together
datalist<-list()
for (i in 1:length(TransFiles)){
  if (i==1|i==2|i==3|i==4){
    mdata<-read.csv(paste0("Output/MutFreq_fil/",TransFiles[i]),stringsAsFactors=FALSE,row.names=1)
    mdata<-mdata[1:NofSamples,1:2]  #use A and T only for nonCpGmaking mutations
    dat<-melt(mdata)
    dat<-dat[!(is.na(dat$value)),]
    
    if (i==1) dat$variable<-c(rep("NonSyn_CpG", n=length(dat$variable)))
    if (i==2) dat$variable<-c(rep("Syn_CpG", n=length(dat$variable)))
    if (i==3) dat$variable<-c(rep("NonSyn_nonCpG", n=length(dat$variable)))
    if (i==4) dat$variable<-c(rep("Syn_nonCpG", n=length(dat$variable)))
  }
  if (i==5|i==6){
    mdata<-read.csv(paste0("Output/MutFreq_fil/",TransFiles[i]),stringsAsFactors=FALSE,row.names=1)
    mdata<-mdata[1:NofSamples,]
    dat<-melt(mdata)
    dat<-dat[!(is.na(dat$value)),]
    if (i==5) dat$variable<-c(rep("NonSyn", n=length(dat$variable)))
    if (i==6) dat$variable<-c(rep("Syn", n=length(dat$variable)))
  }
  datalist[[i]]<-dat
}

Transitions<-do.call(rbind,datalist)        
title<-"Transition Mutations"
MFplot<-ggplot(Transitions,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Mutation Type",y="Mutation frequency")+
  scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC","#4477AACC")) + theme_classic()+
  guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
  theme(axis.text.y = element_text(size =10))+
  ggtitle(paste0(title)) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0("Output/MutFreq_fil/",title,"_Summary(A,T only).pdf"),width=9, height=7, units='in',plot=MFplot)

dev.off()



####
Tv1Files<-Tv1Files[-6]
Tv2Files<-Tv2Files[-6]

datalist.tv<-list()
for (i in 1:length(Tv1Files)){
  if (i==1|i==2){
    mdata1<-read.csv(paste0("Output/MutFreq_fil/",Tv1Files[i]),stringsAsFactors=FALSE,row.names=1)
    mdata2<-read.csv(paste0("Output/MutFreq_fil/",Tv2Files[i]),stringsAsFactors=FALSE,row.names=1)
    mdata1<-mdata1[1:NofSamples,]
    mdata2<-mdata2[1:NofSamples,]
    mdata<-cbind(mdata1[,c(1,4)],mdata2[,c(2,3)])
    dat<-melt(mdata)
    if (i==1) dat$variable<-c(rep("NonSyn_CpG", n=length(dat$variable)))
    if (i==2) dat$variable<-c(rep("Syn_CpG", n=length(dat$variable)))
  }
  
  else {
    mdata1<-read.csv(paste0("Output/MutFreq_fil/",Tv1Files[i]),stringsAsFactors=FALSE,row.names=1)
    mdata2<-read.csv(paste0("Output/MutFreq_fil/",Tv2Files[i]),stringsAsFactors=FALSE,row.names=1)
    mdata1<-mdata1[1:NofSamples,]
    mdata2<-mdata2[1:NofSamples,]
    mdata<-mdata1+mdata2
    dat<-melt(mdata)
    if (i==3) dat$variable<-c(rep("NonSyn_nonCpG", n=length(dat$variable)))
    if (i==4) dat$variable<-c(rep("Syn_nonCpG", n=length(dat$variable)))
    if (i==5) dat$variable<-c(rep("NonSyn", n=length(dat$variable)))
    if (i==6) dat$variable<-c(rep("Syn", n=length(dat$variable)))
  }
  datalist.tv[[i]]<-dat
}

Transvs<-do.call(rbind,datalist.tv)        
title<-"Transversion Mutations"
MFplot2<-ggplot(Transvs,aes(x=variable,y=value,fill=variable))+geom_boxplot(alpha=0.5)+labs(x="Mutation Type",y="Mutation frequency")+
  scale_fill_manual(values=c("#66CCEECC","#228833CC","#CCBB44CC","#EE6677CC","#AA3377CC","#4477AACC")) + theme_classic()+
  guides(fill=FALSE) +theme(axis.text.x = element_text(size =10))+
  theme(axis.text.y = element_text(size =10))+
  ggtitle(paste0(title)) + theme(plot.title = element_text(lineheight=3, face="bold", color="black", size=13))+
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename=paste0("Output/MutFreq_fil/",title,"_Summary.pdf"),width=9, height=7, units='in',plot=MFplot2)

dev.off()
