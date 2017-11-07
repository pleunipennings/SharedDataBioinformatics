##Deshawn's code for syn/non-syn/nonsense

library(seqinr)
library(ape)

#Read in shortened entero data
DF<-read.fasta("InfluenzaBvirus_NA.fasta.mu.trim05")

#Create the data frame
Pos<-c(1:1401)
DF<-data.frame(Pos)
DF$wtnt=""
DF$TrNtFreq=""
DF$WTAA=""
DF$MUTAA=""
DF$WTAAcat=""
DF$MUTAAcat=""
DF$bigAAchange=""
DF$TypeOfSite=""


#Insert Synonymous/Nonsynonymous data into dataframe
####Determining the category of S/N (TypeOfSite)
 #Assigning Synonymous/Nonsynonymous Function
 
   for(i in 1:length(DF$WTAA)){
              if(DF$WTAA[i] == DF$MUTAA[i]){DF$TypeOfSite[i] <- "syn"}
              if(DF$WTAA[i] != DF$MUTAA[i]){DF$TypeOfSite[i] <- "non"}
            # DF$TypeOfSite<- length(which(DF$WTAA==DF$MUTAA)){return(syn)}
           # DF$TypeOfSite<- length(which(DF$WTAA=!DF$MUTAA)){return(non)}
      }
#DF$TypeOfSite<-Syn(s)