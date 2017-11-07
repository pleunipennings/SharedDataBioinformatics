###AMIR JABERI From 11 AM Section 
###OTHER GROUP MEMBERS:Annie Shieh, Nathan O'niel, Andrew Mahoney 
###Many of the coding was done by Nathan O'niel, but I have tried to annotate in my own words as much as possible. 

###This function is for figure two of manuscript###
###This function will allow us to plot the frequency of CpG and non CpG sequences in a given data frame###
###This function will also plot the location of the CpG sequences###


#set working directory
setwd("~/Desktop/GITHUB_FOLDER/SharedDataBioinformatics/11AM_Influ_team")


#Load packages
#Create data vidualizations using gramme of graphics 
library(ggplot2)

#Tools for splitting, applying and combining data 
library(plyr)

#Allows us to have grids graphics in plot
library(grid)

#scale functions for visualization 
library(scales)

#"miscellaneous functions for "grid" graphics 
library(gridExtra)

#Biological sequences retreval and analysis 
library(seqinr)


#Read alignment from the designated file into HPIV1a
HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")

#HPIVIa is a an arbitrary name of a csv file that was made containing the sequence data
HPIV1a

#Gives the consensus sequences for data in file HPIV1a
seqinr::consensus(HPIV1a)

#Gives you length of the consensus sequence. 
length(seqinr::consensus(HPIV1a))

#We know what the length of consensus sequence is, therefore we can create... 
#...position from 1 to the determined length of sequence. 
Pos<-c(1:15473)

#Create data frame called HPIVdf with following variables: 
#The variables will be empty for now. 
HPIVdf<-data.frame(Pos)
HPIVdf$WtNt=""
HPIVdf$MeanFreq=""
HPIVdf$WTAA=""
HPIVdf$MUTAA=""
HPIVdf$WTAAcat=""
HPIVdf$MUTAAcat=""
HPIVdf$bigAAchange=""
HPIVdf$MakesCpG=""
HPIVdf$TypeOfSite=""

#Check to see how data frame looks. 
head(HPIVdf)
tail(HPIVdf)

#place consensus sequence into WtNt variable of data frame HPIVdf
HPIVdf$WtNt<- seqinr::consensus(HPIV1a)

# value "n" will represent our data.frame of use
CpG_Plotter <- function(n){ 
  
  #this adds one value in "n" data.frame to the "0" and "1" place holder values into "1" and "2".
  #this is for the as.integer condition to work
  if (T) {n$makesCpG <- n$makesCpG+3}
  
  #lists which variables return a "2" these make a CpG island when mutated "yes cpg or Y" 
  YCpG <- which(n$makesCpG=="4") 
  
  #lists which variables return a "1" these do not make a CpG island when mutated "No cpg or N"
  NCpG <- which(n$makesCpG=="3") 
  
  #create a Value that looks at mean frequency by yes cpg  
  x1 <- n$MeanFreq[YCpG] 
  
  #create a Value that looks at mean frequency by no cpg
  x2 <- n$MeanFreq[NCpG] 
  
  #create a Value that looks at mean frequency by no cpg
  
  plot.default( x = c(x1, x2),
                xlab = "Position of nucleotide", ylab = "Frequency", 
                main = "Position vs Frequency CpG non-CpG Plot",
                
                #colors the values by their integer, using default association.
                col = (as.integer(n$makesCpG)),  
                
                #the plot is put into a log scale 
                log = "y",)
  
  legend("topright", legend=c("Non-CpG", "CpG"),
         fill =c("green","blue"), cex=0.9,
         text.font=3, bg='lightgray' )
  
  
  return()
}






