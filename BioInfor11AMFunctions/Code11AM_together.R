
#Datafiles<-list.files("OriginalData/") #which data do we have

#Get all data files 
Datafiles<-list.files(recursive = TRUE, full.names = TRUE, pattern = "trim")

filename = Datafiles[1] #start with the first file

print(filename)

source("BioInfor11AMFunctions/meanFreq.R")
DF<-meanFreq(filename)

#look at the protein
paste(translate(DF$wtnt), collapse = "")
#Note: "./11AM_Influ_team/BKpolyomavirus_VP1.fasta.mu.trim05.txt" has the VP1 gene. It is in frame :-) 

source("BioInfor11AMFunctions/getWTAA.R")
DF<-getWTAA(DF)

#look at the protein in WTAA
paste(DF$WTAA[seq(1,nrow(DF),by=3)], collapse = "")

source("BioInfor11AMFunctions/getMUTAA.R")
DF<-getMUTAA(DF)

source("BioInfor11AMFunctions/big_aa_change.R") # Hey, this seems to be case insensitive! I had no idea 
DF<-big_aa_change(DF)

source("BioInfor11AMFunctions/functionSynNonSyn.R")
DF<-functionSynNonSyn(DF)

source("BioInfor11AMFunctions/CpG_finder.R")
DF<-CpG_finder(DF)

DF2<-DF

source("BioInfor11AMFunctions/plotCpGByFeq&Location.R")
LvsF_CpG_Printer(DF)

source("BioInfor11AMFunctions/figure2_plot_meanfreq_syn_nonsyn.R")
plotsyn(DF,"Title",logy)

source("BioInfor11AMFunctions/figureThree.R")
Fig3(DF,"Bocavirus")

source("9AMBioinformaticsFunctions/Figure3_9AM.R")
comparing_mutation_graph(DF)
DF$bigAAchange<-as.factor(DF$bigAAchange)
