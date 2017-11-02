#Script to analyse the frequency data and associate with features. 
#Using the Bacheler et al data

#This script does / makes
#Tests whether non syn muts, syn muts and nonsense muts are different in freq
#Make SingleSiteFrequencySpectraPRO_58.pdf

#* Read the csv files 
#* Perform necessary calcuations
#* Plot results (eventually new script)

#Prep data
if (TRUE){
  library(scales)
  library(plotrix)
  library(RColorBrewer)
  library(sfsmisc)
  #July 2017 now read freqPatTs_Bacheler_Threshold05.csv  and OverviewSelCoeff_BachelerFilter.csv
  #read.table("../Output/freqPatTs_Bacheler_Threshold1.csv",sep=",",header=TRUE,row.names=1)->freqPatTs0
  read.csv("OverviewSelCoeff_BachelerFilter.csv")->OverviewDF
}
#Test whether non syn muts, syn muts and nonsense muts are different in freq

FreqsSyn<-OverviewDF$MeanFreq[OverviewDF$TypeOfSite=="syn"]
FreqsNonSyn<-OverviewDF$MeanFreq[OverviewDF$TypeOfSite=="nonsyn"]
FreqsStop<-OverviewDF$MeanFreq[OverviewDF$TypeOfSite=="stop"]

wilcox.test(FreqsSyn, FreqsNonSyn,alternative = "greater", paired = FALSE)
wilcox.test(FreqsNonSyn,FreqsStop,alternative = "greater", paired = FALSE)

# ?wilcox.test
#Make a figure with the selection coefficients across Pol
#Currently figure 2 in the revision P Genetics Sept 2017
if (TRUE){
  #pdf("../Output/EstSelCoeffPRO_aug2017.pdf",width=12,height=8)
  png("EstSelCoeffPRO_aug2017.png",width=12,height=7.5,units="in",res=100)
#?png #prints out the below plot as a png!
  par(mfrow=c(1,1))
  maxnuc=1000
#??maxnuc #maxnuc? maxnuc sets the maximum number of variables to be used in 
  # the plot
  par(mar = c(3,5,1,2))
#?par #what does par do?
  plot(OverviewDF$num[40:maxnuc],OverviewDF$EstSelCoeff[40:maxnuc],
       log="y", ylab="Estimated Selection Coefficient (cost)",cex.lab=1.3,
       xaxt="n",yaxt="n", xlab="",
       col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-4,1),xlim=c(40,979))
  
  axis(1,at=c(3*seq(15,95,by=20)-1,296+30),labels=c(seq(15,95,by=20),""))
  axis(1,at=3*seq(109,349,by=20)-1,labels=seq(109-99,349-99,by=20))
  axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),las=1,line=0,tick=FALSE)
  eaxis(side = 2, at = 10^((-0):(-(5))),label=rep("",6))
  
  #color Protease region grey
  rect(0, 0.00001, 297.5, 2, density = NULL, angle = 45,col="grey70",border = NA)
  for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray41")}
  
  cols <- brewer.pal(6, "Set2")[c(1, 2, 3, 6)]
  for (i in 40:maxnuc){
    c=0; co = 1
    if (OverviewDF$TypeOfSite[i]=="stop"&OverviewDF$WTnt[i]%in%c("g","c")) {c=1;p=21}
    if (OverviewDF$TypeOfSite[i]=="syn"&OverviewDF$WTnt[i]%in%c("g","c")) {c=cols[1];p=21}
    if (OverviewDF$TypeOfSite[i]=="syn"&OverviewDF$WTnt[i]%in%c("a","t")) {c=cols[1];p=21}
    if (OverviewDF$TypeOfSite[i]=="nonsyn"&OverviewDF$WTnt[i]%in%c("c","g")) {c=cols[2];p=21}
    if (OverviewDF$TypeOfSite[i]=="nonsyn"&OverviewDF$WTnt[i]%in%c("a","t")) {c=cols[4];p=21}
    if (i %in% 73:81) {p = 22; co = 2} #for Active site Protease change pch
    if (c!=0) points(OverviewDF$num[i],OverviewDF$EstSelCoeff[i],pch=p,col=co,
                     bg=rgb(red=col2rgb(c)[1]/255,
                            green=col2rgb(c)[2]/255,
                            blue=col2rgb(c)[3]/255,
                            maxColorValue = 1,alpha=0.8),
                     cex=2)
  }
  
  #Add "Protease" and "RT" words
  rect(0, 0.000001, 1200, 3.5*10^-4, density = NULL, angle = 45,col=1,border = NA)
  text(55*3,2.9*10^-4,"PROTEASE",col="white")
  rect(297.5, 0.000001, 1200, 3.5*10^-4, density = NULL, angle = 45,col="grey40",border = NA)
  text(220*3,2.9*10^-4,"REVERSE TRANSCRIPTASE",col="white")
  
  
  
  #Add legend
  legpos=296;legposV=0.4
  rect(legpos*3, 0.4*legposV, (legpos+42.5)*3, 1.7*legposV, density = NULL, angle = 45,col=alpha("white",1))
  points((legpos+5)*3,legposV/0.7,pch=21,bg=1,col=1,cex=2)
  text((legpos+9)*3,legposV/0.7,"Nonsense",adj=0)
  points((legpos+5)*3,legposV,pch=21,bg=cols[2],col=1,cex=2)
  text((legpos+9)*3,legposV,"Non-syn, C/G",adj=0)
  points((legpos+5)*3,legposV*0.7,pch=21,bg=cols[4],col=1,cex=2)
  text((legpos+9)*3,legposV*0.7,"Non-syn, A/T",adj=0)
  points((legpos+5)*3,legposV*0.49,pch=21,bg=cols[1],col=1,cex=2)
  text((legpos+9)*3,legposV*0.49,"Synonymous",adj=0)
  
  dev.off()
}

#####Make a figure with single site frequency spectra for Protease AA 58####
#Currently figure 1 in the revision P Genetics Sept 2017
if (FALSE){
  pdf("../Output/SingleSiteFrequencySpectraPRO_58_July2017.pdf",width=8,height=4)
  zerobar=50; h2=22; x1=0.25
  cols <- c(0,brewer.pal(6, "Set2")[c(2, 1)])
  par(mfrow=c(2,3))
  par(mar = c(1,3,4,2))
  for (i in 172:174){
    #first create empty plot with title
    if (i == 172){
      #par(fig=c(0,2/3,0,1))
      t=paste("nonsense mutation",sep="")
      hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
           col=cols[1],border=0,
           #    main= bquote(paste(.(t),(C %->% T ))), cex=1.3,
           main="",cex=1.2,
           xlab="Frequency", ylab="Count",cex.lab=1.4)
      title(t,cex=1.2,line=0)
      text(x1,30,"observed data",cex=1.3)
      text(x1,h2,"(C172T)",cex=1.3)
    }
    if (i == 173){
      t=paste("         non-synonymous mutation",sep="")
      #    t=paste("Protease: site ", i,"\n non-synonymous mutation",sep="")
      hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
           col=cols[2],border=0,
           #    main = bquote(paste(.(t),(A %->% G ))), cex=1.3,
           main= "", cex=1.2,
           xlab="Frequency", ylab="Count",cex.lab=1.4)
      title(t,cex=1.2,line=0)
      text(x1,30,"observed data",cex=1.3)
      text(x1,h2,"(A173G)",cex=1.3)
      
    }
    if (i == 174){
      t=paste("   synonymous mutation",sep="")
      #    t=paste("Protease: site ", i,"\n synonymous mutation",sep="")
      hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
           col=cols[3],border=0,
           #    main = bquote(paste(.(t),(G %->% A ))), cex=1.3,
           main= "", cex=1.2,
           xlab="Frequency", ylab="Count",cex.lab=1.4)
      title(t,cex=1.2,line=0)
      text(x1,30,"observed data",cex=1.3)
      text(x1,h2,"(G174A)",cex=1.3)
      
    }
    #Next, show  0 bar
    if (i == 172){
      hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
           yaxt="n",col=OverviewDF$color[which(OverviewDF$num==i)],add=T)}
    if (i == 173){
      hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
           yaxt="n",col=cols[2],add=T)}
    if (i == 174){
      hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
           yaxt="n",col=cols[3],add=T)}
    
    #next show all data (unfiltered), but only until 50 for 0 cat
    if (i == 172){
      hist(c(rep(0,min(zerobar-10,length(which(freqPatTs0[,i]<0.02)))),freqPatTs0[,i][which(freqPatTs0[,i]>0)]),
           breaks=seq(0,1,by=0.02),add=T,
           col=OverviewDF$color[which(OverviewDF$num==i)])}
    if (i == 173){
      hist(c(rep(0,min(zerobar-10,length(which(freqPatTs0[,i]<0.02)))),freqPatTs0[,i][which(freqPatTs0[,i]>0)]),
           breaks=seq(0,1,by=0.02),add=T,
           col=cols[2])}
    if (i == 174){
      hist(c(rep(0,min(zerobar-10,length(which(freqPatTs0[,i]<0.02)))),freqPatTs0[,i][which(freqPatTs0[,i]>0)]),
           breaks=seq(0,1,by=0.02),add=T,
           col=cols[3])}
    
    axis(2,labels = c(10,20,30,max(zerobar,length(which(freqPatTs0[,i]<0.02)))), 
         at = c(10,20,30,zerobar), las=1)
    if (length(which(freqPatTs0[,i]<0.02))>=zerobar){
      axis.break(axis=2,breakpos=zerobar-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
      points(c(0.01,0.02),c(zerobar-10,zerobar-10),pch=15,cex=2.5,col="white")
    }else{axis(2,labels = zerobar-10,at=zerobar-10,las=1)}
    
  }
  #next, show simulated frequencies
  #    par(fig=c(2/3,1,0,1), new = TRUE)
  par(mar = c(4.5,3,1,2))
  for (i in 172:174){
    if (i ==172)Freqs<-read.csv("../Output/SimFreqs172.csv",row.names=1)[1][,1]
    if (i ==173)Freqs<-read.csv("../Output/SimFreqs173.csv",row.names=1)[1][,1]
    if (i ==174)Freqs<-read.csv("../Output/SimFreqs174.csv",row.names=1)[1][,1]
    t=paste("simulated data",sep="")
    hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),yaxt="n",
         col=cols[1],border=0,
         main="",cex=1.2,
         xlab="Frequency", ylab="Count",cex.lab=1.4)
    #title(t,cex=1.2,line=0)
    text(x1,30,"simulated data",cex=1.3)
    if (i ==172)text(x1,h2,"(s=1)",cex=1.3)
    if (i ==173)text(x1,h2,paste("(s=",round(OverviewDF$EstSelCoeff[173],3),")",sep=""),cex=1.3)
    if (i ==174)text(x1,h2,paste("(s=",round(OverviewDF$EstSelCoeff[174],3),")",sep=""),cex=1.3)
    
    hist(rep(0,zerobar),breaks=seq(0,1,by=0.02),xlim=c(0,.5),ylim=c(0,zerobar),
         yaxt="n",col=brewer.pal(4, "Set2")[3],add=T)
    hist(c(rep(0,
               min(zerobar,length(which(Freqs<0.02)))
    ),
    Freqs[which(Freqs>=0.02)]),
    breaks=seq(0,1,by=0.02),add=T,
    col=brewer.pal(4, "Set2")[3])
    axis(2,labels = c(10,20,30,max(zerobar,length(which(Freqs<0.02)))), 
         at = c(10,20,30,zerobar), las=1)
    if (length(which(Freqs<0.02))>=zerobar){
      axis.break(axis=2,breakpos=zerobar-10,bgcol="white",breakcol="black",style="slash",brw=0.02)
      points(c(0.01,0.02),c(zerobar-10,zerobar-10),pch=15,cex=2.5,col="white")
    }else{axis(2,labels = zerobar-10,at=zerobar-10,las=1)}
  }
  
  dev.off()
}