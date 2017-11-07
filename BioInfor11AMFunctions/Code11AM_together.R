

source("BioInfor11AMFunctions/meanFreq.R")
DF<-meanFreq("denguesmall.txt")
DF2<-DF
DF<-DF2

source("BioInfor11AMFunctions/getWTAA.R")
DF<-getWTAA(DF)

source("BioInfor11AMFunctions/getMUTAA.R")
DF<-getMUTAA(DF)

source("BioInfor11AMFunctions/big_aa_change.R")
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
Fig3(DF)
