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
levels(data$number) <- gsub("a", as.numeric("1"), levels(data$number))
levels(data$number) <- gsub("c", as.numeric("2"), levels(data$number))
levels(data$number) <- gsub("g", as.numeric("3"), levels(data$number))
levels(data$number) <- gsub("t", as.numeric("4"), levels(data$number))

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

syn <- which(data$TypeOfSite=="syn")
non <- which(data$TypeOfSite == "nonsyn")
cols<-c("red","yellow","blue","green")
colsyn<-cols[syndata$combo]


#syn subset stuff
syndata <- subset(data, TypeOfSite=="syn")
noncpGdata <- subset(syndata, combo=="noAA noCPG")
cpGdata<- subset(syndata, combo=="noAA yesCPG")


#building the subset plot
ggplot(aes(factor(WTnt), MeanFreq), data = noncpGdata)+
    geom_jitter(col = "red") 
geom_jitter(aes(factor(WTnt), MeanFreq), data = cpGdata)+
    geom_jitter(color ="blue")





#nonsyn sub set stuff
nonsyndata <- subset(data, TypeOfSite=="nonsyn")
nonnoncpGdata <- subset(nonsyndata, combo=="noAA noCPG")
noncpGdata<- subset(nonsyndata, combo=="noAA yesCPG")



# not usinf ggplot
dotchart( as.numeric(syndata$num),syndata$MeanFreq, type="p", 
          prob = TRUE, col=col)


plot(jitter(as.numeric(syndata$WTnt)), syndata$MeanFreq, type="p", prob = TRUE, col=colsyn, pch=16)
boxplot(MeanFreq~as.numeric(WTnt), data=syndata, col=colsyn)

?boxplot()

#hist(pdg$Course.total[IndsMen], prob = TRUE, breaks = 30, col= rgb(0,0.,1,0.5), add = TRUE)


#graph 2 Non-synomymous Sites
