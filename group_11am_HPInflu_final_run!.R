#######
HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
HPIVdf<-data.frame(Pos)
HPIVdf$WtNt=""
HPIVdf$WtNt<- seqinr::consensus(HPIV1a)
head(HPIVdf)

df <- HPIVdf
df

source("BioInfor11AMFunctions/meanFreq.R")
df <- meanFreq(df)

source("BioInfor11AMFunctions/getWTAA.R")
df <- getWTAA(df)

source("BioInfor11AMFunctions/getMUTAA.R")
df <- getMUTAA(df)

source("BioInfor11AMFunctions/Big_AA_Change.R")
df <- Big_AA_Change(df)



source("BioInfor11AMFunctions/functionSynNonSyn.R")
df <- functionSynNonSyn(df)

source("BioInfor11AMFunctions/CpG_finder.R")
CpG_finder <- CpG_finder(df)
head(df)
#plots

source("BioInfor11AMFunctions/plotCpGByFeq&Location.R")



######################################
source("BioInfor11AMFunctions/meanFreq.R")
df<-meanFreq("11AM_Influ_team/humanparainfluenzavirus3_HN.fasta_pruned.mu.trim05")
df2<-df
df<-df2
head(df)

source("BioInfor11AMFunctions/functionWTAA.R") #stuck here
df$WTAA=""
df<-getWTAA(df)
df <- df[1:(nrow(df)-2),]
head(df)
source("BioInfor11AMFunctions/getMUTAA.R")
df<-getMUTAA(df)

source("BioInfor11AMFunctions/big_aa_change.R")
df<-big_aa_change(df)

source("BioInfor11AMFunctions/functionSynNonSyn.R")
df<-functionSynNonSyn(df)

source("BioInfor11AMFunctions/CpG_finder.R")
df<-CpG_finder(df)

df2<-df
head(df)
source("BioInfor11AMFunctions/plotCpGByFeq&Location.R")
LvsF_CpG_Printer(df)

source("BioInfor11AMFunctions/figure2_plot_meanfreq_syn_nonsyn.R")
plotsyn(df,"Title",logy)

source("BioInfor11AMFunctions/figureThree.R")
Fig3(df)
