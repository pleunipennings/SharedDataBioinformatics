# testerino



#group 6 Victoria, Sam, Mordy, Gordan
#Metting Wednesday 7-9 
# BK polyom/na
#Due October 31

setwd("~/Desktop/Git")

install.packages('ggplot2')
install.packages("gridExtra")
install.packages("dplyr")
library(ggplot2)
library(gridExtra)
library(dplyr)



data<-read.csv('OverviewSelCoeff_BachelerFilter.csv') 

#Example
library(ggplot2)

?theme_set()
# test plot sample
#g <- ggplot(mpg, aes(manufacturer, cty))
#g + geom_boxplot() + 
#    geom_dotplot(binaxis='y', 
#                 stackdir='center', 
#                 dotsize = .5, 
 #                fill="red") +
 #   theme(axis.text.x = element_text(angle=65, vjust=0.6)) + 
#    labs(title="Box plot + Dot plot", 
#         subtitle="City Mileage vs Class: Each dot represents 1 row in source data",
#         caption="Source: mpg",
#         x="Class of Vehicle",
#         y="City Mileage")


#building the combo lines to help sort 

data$combo<- (data$bigAAChange*3) + (data$makesCpG*2)
data$combo<- as.factor(data$combo)
levels(data$combo) <- gsub("0", "noAA noCPG", levels(data$combo))
levels(data$combo) <- gsub("2", "noAA yesCPG", levels(data$combo))
levels(data$combo) <- gsub("3", "yesAA noCPG", levels(data$combo))
levels(data$combo) <- gsub("5", "yesAA yesCPG", levels(data$combo))

#building numbers for letter might not be needed
data$number<- data$WTnt
data$number<- as.factor(data$num)

levels(data$number) <- gsub("a", as.numeric("1"), levels(data$num))
levels(data$number) <- gsub("c", as.numeric("2"), levels(data$num))
levels(data$number) <- gsub("g", as.numeric("3"), levels(data$num))
levels(data$number) <- gsub("t", as.numeric("4"), levels(data$num))

syn <- which(data$TypeOfSite=="syn")
non <- which(data$TypeOfSite == "nonsyn")
cols<-c("red","yellow","blue","green")
#colsyn<-cols[syndata$combo]

#syn subset stuff
syndata <- subset(data, TypeOfSite=="syn")

#syn subset for no AA and no CPG
synNNdata <- subset(syndata, combo=="noAA noCPG")
#syn subset for no AA and no CPG for a, c, g, t (for HIV all 4 should be here)
synNNa <- subset(synNNdata, WTnt=="a")
synNNc <- subset(synNNdata, WTnt=="c")
synNNg  <- subset(synNNdata, WTnt=="g")
synNNt  <- subset(synNNdata, WTnt=="t")
#syn subset for no AA and yes CPG
synNYdata<- subset(syndata, combo=="noAA yesCPG")
#syn subset for no AA and yes CPG for a, c, g, t (for HIV a, t )
synNYa <- subset(synNYdata, WTnt=="a")
synNYc <- subset(synNYdata, WTnt=="c")
synNYg  <- subset(synNYdata, WTnt=="g")
synNYt  <- subset(synNYdata, WTnt=="t")
#syn subset for yes AA and no CPG 
synYNdata <- subset(syndata, combo=="yesAA noCPG")
#syn subset for yes AA and no CPG for a, c, g, t (none for HIV)
synNNa <- subset(synYNdata, WTnt=="a")
synNNc <- subset(synYNdata, WTnt=="c")
synNNg  <- subset(synYNdata, WTnt=="g")
synNNt  <- subset(synYNdata, WTnt=="t")
#syn subset for yes AA and yes CPG
synYYdata <- subset(syndata, combo=="yesAA yesCPG")
#syn subset for yes AA and yes CPG for a, c, g, t (none for HIV)
synYYa <- subset(synYYdata, WTnt=="a")
synYYc <- subset(synYYdata, WTnt=="c")
synYYg  <- subset(synYYdata, WTnt=="g")
synYYt  <- subset(synYYdata, WTnt=="t")


#graph 1 Synonymous Sites

#still working
graph <- ggplot(aes(factor(WTnt), MeanFreq, colour = combo), data = syndata)+
    geom_jitter()
#not yet done
graph + geom_boxplot() + 
    geom_dotplot(binaxis='y', 
                 stackdir='center', 
                 dotsize = .5, 
                 fill="blue") +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 






fora <- which(noncpGdata$WTnt=="a")

#building the subset plot
ggplot(aes(factor(WTnt), MeanFreq), data = noncpGdata)+
    geom_jitter(col = "red") +
    geom_errorbar(aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(aes('a',median(c(median(lowerConf),median(upperConf)))))
# need to find lower and upper conf for each individual nuc (a,c,g,t)

ggplot(aes(factor(WTnt), MeanFreq), data = cpGdata)+
    geom_jitter(color ="blue")





#nonsyn sub set stuff
nonsyndata <- subset(data, TypeOfSite=="nonsyn")

nonyescpGdata<- subset(nonsyndata, combo=="noAA yesCPG")

#nonsyn subset for no AA and no CPG
nonNNdata <- subset(nonsyndata, combo=="noAA noCPG")
#nonsyn subset for no AA and no CPG for a, c, g, t (for HIV all 4 should be here)
nonsynNNa <- subset(nonNNdata, WTnt=="a")
nonsynNNc <- subset(nonNNdata, WTnt=="c")
nonsynNNg  <- subset(nonNNdata, WTnt=="g")
nonsynNNt  <- subset(nonNNdata, WTnt=="t")
#nonsyn subset for no AA and yes CPG
nonNYdata<- subset(nonsyndata, combo=="noAA yesCPG")
#syn subset for no AA and yes CPG for a, c, g, t (for HIV a, t )
nonsynNYa <- subset(nonNYdata, WTnt=="a")
nonsynNYc <- subset(nonNYdata, WTnt=="c")
nonsynNYg  <- subset(nonNYdata, WTnt=="g")
nonsynNYt  <- subset(nonNYdata, WTnt=="t")
#syn subset for yes AA and no CPG 
nonYNdata <- subset(nonsyndata, combo=="yesAA noCPG")
#syn subset for yes AA and no CPG for a, c, g, t (for HIV all 4 should be here)
nonsynNNa <- subset(nonYNdata, WTnt=="a")
nonsynNNc <- subset(nonYNdata, WTnt=="c")
nonsynNNg  <- subset(nonYNdata, WTnt=="g")
nonsynNNt  <- subset(nonYNdata, WTnt=="t")
#syn subset for yes AA and yes CPG
nonYYdata <- subset(nonsyndata, combo=="yesAA yesCPG")
#syn subset for yes AA and yes CPG for a, c, g, t (for HIV a, t)
nonsynYYa <- subset(nonYYdata, WTnt=="a")
nonsynYYc <- subset(nonYYdata, WTnt=="c")
nonsynYYg  <- subset(nonYYdata, WTnt=="g")
nonsynYYt  <- subset(nonYYdata, WTnt=="t")


# not usinf ggplot
dotchart( as.numeric(syndata$num),syndata$MeanFreq, type="p", 
          prob = TRUE, col=col)


plot(jitter(as.numeric(syndata$WTnt)), syndata$MeanFreq, type="p", prob = TRUE, col=colsyn, pch=16)
boxplot(MeanFreq~as.numeric(WTnt), data=syndata, col=colsyn)

?boxplot()

#hist(pdg$Course.total[IndsMen], prob = TRUE, breaks = 30, col= rgb(0,0.,1,0.5), add = TRUE)


#graph 2 Non-synomymous Sites
