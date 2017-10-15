# Base Scripts
#Working on this Sep 2015, preparing code for Marion and Kristof

#Working on this again in December 2014. Want to test how well we can estimate the mean frequency of a mutation. 
#I will determine mean freq for each site and then do bootstrapping to get 95% conf intervals
#This may be useful as prelim data for NSF proposal. 

#I plan to plot the frequencies of mutations WITHIN all patients. 
#So for each site, I determine the B consensus base, this will be WT
#Next, for each patient, I determine whether the seqs on day 1 were WT. 
#If that is the case then I will look at the freq of the non-WT bases on all days after day 1. 
#Accross all patients, I will have around one hundred frequencies for each base and each site. Now I look at the freq dist for each site. 
#If all 4-fold sites are truly neutral, then they should all show the same frequency distribution. 
#The distribution will not look like a neutral one because of the conditioning on starting off entirely WT. But it should still be OK to compare between sites, because I do the conditioning for each site.

#DEPENDS ON source("/Users/pleuni/Documents/Research/HIV/SoftSweepsInHIV/Bacheler2000/RResistanceMutations.r")
#DEPENDS ON "/Users/pleuni/Documents/Research/HIV/SoftSweepsInHIV/HowToGenbank/HIV1_CON_2004_POL_DNA.fasta"
#DEPENDS ON OR CREATES "freqPatSite.csv" (No longer, now this is made in createFrequencies-Bacheler.R)
#CREATES pdf("Distribution_Prob_Seg.pdf")
#CREATES pdf("freqdis_WITHINallpatients_12.pdf")
#CREATES pdf("freqdis_WITHINallpatients_3_ffdeg.pdf")
#CREATES pdf("freqdis_WITHINallpatients_3_NOTffdeg.pdf")

#load relevant libraries and read consensusfasta file
CurrentWD<-getwd()
setwd("~/Documents/Git/bachelerProject/Rscripts")

library(ape)
library(seqinr)
library(pegas)
#read the file with the resistance mutations
source("./RResistanceMutations.r")
#read the fasta file 
consensusfasta<-read.dna("../Data/HIV1_CON_2004_POL_DNA.fasta", format = "fasta",as.character=TRUE)	
#where is the start of POL? 
polstart=regexpr("cctca",paste(consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_B"),],collapse=""))[1]
consensusB<-consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_B"), polstart:(polstart+983)]
consensusA<-consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_A1"), polstart:(polstart+983)]
consensusC<-consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_C"), polstart:(polstart+983)]
consofcons<-consensusfasta[which(row.names(consensusfasta)=="CON_OF_CONS"), polstart:(polstart+983)]
consensus01AE<-consensusfasta[which(row.names(consensusfasta)=="CONSENSUS_01_AE"), polstart:(polstart+983)]
list.files(path="../Data/BachelerFiles/FASTAfiles/")->listfastafiles

setwd(CurrentWD)

#* Transition function*
transition<-function(nuc){
  if (nuc=="a") return("g")
  if (nuc=="g") return("a")
  if (nuc=="c") return("t")
  if (nuc=="t") return("c")
}

typeofsitefunction<-function(WTcodon, mutantcodon){
  WTAA<-seqinr::translate(WTcodon)
  MUTAA<-seqinr::translate(mutantcodon)
  if (WTAA == MUTAA) return ("syn")
  else if (MUTAA == "*") return ("stop")
  else return ("nonsyn")
}

TypeOfSite<-c()
for (codon in 1:13) TypeOfSite<-c(TypeOfSite,c("overlap","overlap","overlap"))
for (codon in 14:(984/3)){#for each codon in the sequence
  positions <- c(codon*3-2,codon*3-1, codon*3)
  WTcodon <- consensusB[positions]
  mutant1codon <- c(transition(WTcodon[1]), WTcodon[2:3])
  mutant2codon <- c(WTcodon[1],transition(WTcodon[2]), WTcodon[3])
  mutant3codon <- c(WTcodon[1:2], transition(WTcodon[3]))
  TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant1codon))
  TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant2codon))
  TypeOfSite<-c(TypeOfSite,typeofsitefunction(WTcodon,mutant3codon))
}
#make sure that resistance sites in RT have a diff type of site
TypeOfSite[sort(c((RTImuts$pos*3)-2,(RTImuts$pos*3)-1,(RTImuts$pos*3)))+297]<-"res"
TypeOfSite[sort(c((Indinavirmuts$pos*3)-2,(Indinavirmuts$pos*3)-1,(Indinavirmuts$pos*3)))]<-"res"

EstimatedS <- function(mu, meanfreq){
  if (meanfreq == 0) return (1)
  else return (min(c(mu/meanfreq,1)))
}

#Amino acid changes
pos <- "R|H|K"
neg <- "D|E"
unc <- "S|T|N|Q"
spe <- "C|U|G|P"
hyd <- "A|I|L|F|M|W|Y|V"
amCat <- function(AA){
  if(regexpr(pos, AA) > 0){ return(0) }
  if(regexpr(neg, AA) > 0){ return(1) }
  if(regexpr(unc, AA) > 0){ return(2) }
  if(regexpr(spe, AA) > 0){ return(3) }
  if(regexpr(hyd, AA) > 0){ return(4) }
  return(5)
}