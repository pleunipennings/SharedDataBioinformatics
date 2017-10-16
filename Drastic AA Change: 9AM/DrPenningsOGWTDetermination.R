#This R script reads the Bacheler et al fasta files and creates a csv file that has frequencies for each patient and each site for the transition mutations. 
# We filter the data before writing the file. 

#new addition
#1. only take patients with a consensus sequence at Day0
#2. only consider changed nucleotides when the two other nucleotides in a codon are not mutated

#get all the names of the fasta files for each patient in Bacheler2000 dataset
#determine all the frequencies and store in freqPatSite data.frame
#only needed if the stored data are not O

#Load libraries and necessary files from the baseRscript.Rmd
setwd("~/Documents/Git/bachelerProject")
source('Rscripts/baseRscript.R')

#Read the correct fastafiles.

listfastafiles<-list.files("/Users/hasansulaeman/Google Drive/Bioinformatics/InfluenzaVirus")
lengthlatersequenceslist<-c();lengthallsequenceslist<-c()

###############   
#First, let's remove patients with too few sequences (less than 5). 
Numseqs<-c()
for (i in 1:length(listfastafiles)){ #for each fastafile
  filename=paste("/Users/hasansulaeman/Google Drive/Bioinformatics/InfluenzaVirus",substr(listfastafiles[i],1,6),".fasta",sep="")
  patfasta<-read.dna(filename, format = "fasta",as.character=TRUE) #read the file
  Numseqs<-c(Numseqs,nrow(patfasta))
}
#Remove from listfastafiles the files / patients with too few sequences
if (length(which(Numseqs<5))>0){listfastafiles<-listfastafiles[-which(Numseqs<5)]}

###############
#make dataframe that keeps track, for each position and patients, what happened according to criteria and which data points are included
#filterSummary66<-data.frame(row.names=substr(listfastafiles,1,6))
#WTthreshold = 0.5 # we need at least 66% wt on day 0 to keep the position and the patient in the freqPatTs_threshold

################
#OK, now we go through the fasta files again to get the frequencies

#Pleuni 05/17 continue here, add WTthreshold = 1
#Pleuni July 2017: now think the threshold should be 0.5

for (WTthreshold in c(1)){ #0 means no threshold
  
  #make dataframe with frequencies for all non-muts for all patients for all sites filtered with the WT threshold.
  freqPatTs_threshold<-data.frame(row.names=substr(listfastafiles,1,6))
  CountData<-data.frame(pat=character(), pos=integer(),WTnt=character(),MutCount=integer(),Totalcount=integer(),stringsAsFactors = F)
  Nonconsensusday0_pat_pos<-c()
  
  for (i in 1:length(listfastafiles)){ #for each fastafile
    #for (i in 1:2){ #for each fastafile
    filename=paste("Data/BachelerFiles/FASTAfiles/",substr(listfastafiles[i],1,6),".fasta",sep="")
    print(filename)
    patfasta<-read.dna(filename, format = "fasta",as.character=TRUE) #read the file
    
    # find which seqs are from first day of sampling
    days=sort(unique(as.numeric(substr(row.names(patfasta),5,7))))
    day0sequences<-which(as.numeric(substr(row.names(patfasta),5,7))==days[1])
    latersequences<-which(as.numeric(substr(row.names(patfasta),5,7))>days[1])
    allsequences=c(day0sequences,latersequences)
    
    for (j in 1:984){#for each site in the sequence
      #print(paste('j',j))
      WT =  consensusB[j] #what is WT at site j?
      transitionnuc = transition(WT) #which nuc is the transition mutation?
      #how many of the patfasta[day0sequences,j] are WT?
      fracWT=length(which(patfasta[day0sequences,j]==WT))/length(patfasta[day0sequences,j])
      #print(fracWT)
      if (fracWT<WTthreshold) {
        freqPatTs_threshold[i,j]<-NA
        #keep track of pat and pos that are excluded
        Nonconsensusday0_pat_pos<-rbind(Nonconsensusday0_pat_pos,c(i,j))
      } else if (fracWT>=WTthreshold){
        #check wether the neigboring sequences are the same
        # in order to check whether j is the first, second or third position of the codon, you can do
        if(j %in% seq(1,982,by=3)) {# first position
          #use only those sequences that are WT at pos 2 and 3 of the triplet
          goodsequences<-which(paste(patfasta[,j+1],patfasta[,j+2]) == paste(consensusB[c(j+1)],consensusB[(j+2)]))
          freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]== transitionnuc))/length(goodsequences)
          CountData[length(CountData[,1])+1,]<-c(substr(listfastafiles[i],1,6),j,WT,length(which(patfasta[goodsequences,j]== transitionnuc)),length(goodsequences))
        }
        if((j %in% seq(2,983,by=3))) {# second position
          #use only those sequences that are WT at pos 1 and 3 of the triplet
          goodsequences<-which(paste(patfasta[,j-1],patfasta[,j+1]) == paste(consensusB[c(j-1)],consensusB[(j+1)]))
          freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]==transitionnuc))/length(goodsequences)
          CountData[length(CountData[,1])+1,]<-c(substr(listfastafiles[i],1,6),j,WT,length(which(patfasta[goodsequences,j]== transitionnuc)),length(goodsequences))
        }
        if((j %in% seq(3,984,by=3))) {# third position
          #use only those sequences that are WT at pos 1 and 2 of the triplet
          goodsequences<-which(paste(patfasta[,j-2],patfasta[,j-1]) == paste(consensusB[c(j-2)],consensusB[(j-1)]))
          freqPatTs_threshold[i,j]=length(which(patfasta[goodsequences,j]==transitionnuc))/length(goodsequences)
          CountData[length(CountData[,1])+1,]<-c(substr(listfastafiles[i],1,6),j,WT,length(which(patfasta[goodsequences,j]== transitionnuc)),length(goodsequences))
        }
        if (length(goodsequences)==0)  freqPatTs_threshold[i,j]<-NA
      }
    }
  }
  if (WTthreshold == 0){write.csv(freqPatTs_threshold,file="Output/freqPatTs_Bacheler_Threshold0.csv");write.csv(CountData,"Output/BachelerCountData_Threshold0.csv")}
  if (WTthreshold == 0.5){write.csv(freqPatTs_threshold,file="Output/freqPatTs_Bacheler_Threshold05.csv");write.csv(CountData,"Output/BachelerCountData_Threshold05.csv");write.csv(Nonconsensusday0_pat_pos,"Output/Nonconsensusday0_pat_pos_Threshold05.csv")}
  if (WTthreshold == 1){write.csv(freqPatTs_threshold,file="Output/freqPatTs_Bacheler_Threshold1.csv");write.csv(CountData,"Output/BachelerCountData_Threshold1.csv")}
  
  print(paste("For threshold",WTthreshold))
  print("how much data discarded because of filtering?")
  print(length(Nonconsensusday0_pat_pos[,1])/(length(freqPatTs_threshold[,1])*length(freqPatTs_threshold[1,])))
}

plot(OverviewDFilter$MeanFreq,table(Nonconsensusday0_pat_pos[,1]))

for (i in OverviewDFilter$num){
  OverviewDFilter$FracFiltered[i]<-length(which(Nonconsensusday0_pat_pos[,2]==i))
}
plot(OverviewDFilter$MeanFreq,OverviewDFilter$FracFiltered) 
