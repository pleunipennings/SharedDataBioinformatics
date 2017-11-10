library(seqinr)
library(ape)

#various set directories for file paths in diff parts of my computer
setwd("~/bioinformatics/bioinformaticsproject")
setwd("~/SharedDataBioinformatics/BioInfor11AMFunctions")

fasta_file<-"HumanBocavirus1_VP1.fasta_pruned.mu.trim05"
dengue_basic <- read.fasta(fasta_file)
number_of_seqs <- length(dengue_basic)
dengue_align <- read.alignment(fasta_file, format = "fasta", forceToLower = T)



source("MeanFreq.R")
setwd("~/bioinformatics/bioinformaticsproject")
meanFreq("HumanBocavirus1_VP1.fasta_pruned.mu.trim05")->DF

setwd("~/SharedDataBioinformatics/BioInfor11AMFunctions")
source("functionWTAA.R")
getWTAA(DF)->DF$WTAA

DF->dfwtaa
source("getMUTAA.R")
getMUTAA(dfwtaa)->dfmutwt

source("Big_AA_Change.R")
big_aa_change(dfmutwt)->dfaachange

source("functionSynNonSyn.R")
functionSynNonSyn(dfaachange)-> dfaachange$TypeOfSite
dfaachange->DF1

source("CpG_finder.R")
CpG_finder(DF1)->dffull

#Cpg plot
source("plotCpGByFeq&Location.R")
LvsF_CpG_Printer(dffull)

#plot syn nonsyn
source("figure2_plot_meanfreq_syn_nonsyn.R")
plotsyn(dffull,"Higher Frequency of Synonymous Mutations", logy)

#Figure 3
source("figureThree.R")
Fig3(dffull, "Bocavirus NS1")
