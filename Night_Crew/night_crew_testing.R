# testerino



#group 6 Victoria, Sam, Mordy, Gordan
#Metting Wednesday 7-9 
# BK polyom/na
#Due October 31

setwd("~/Desktop/Git")
setwd("~/Desktop")
install.packages('ggplot2')
library(ggplot2)


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
data$number<- as.factor(data$number)
levels(data$number) <- gsub("a", as.numeric("1"), levels(data$number))
levels(data$number) <- gsub("c", as.numeric("2"), levels(data$number))
levels(data$number) <- gsub("g", as.numeric("3"), levels(data$number))
levels(data$number) <- gsub("t", as.numeric("4"), levels(data$number))

#syn <- which(data$TypeOfSite=="syn")
#non <- which(data$TypeOfSite == "nonsyn")
cols<-c("red","yellow","blue","green")
colsyn<-cols[syndata$combo]

#adding x values
#data$one<- 1
#data$two<- 2
#data$three<- 3
#data$four<- 4
#data$five<- 5
#data$six<-6
#data$seven<- 7
#data$eight<- 8
#data$nine<- 9
#data$ten<- 10
#data$eliven<- 11
#data$twelve<- 12
#data$thirteen<- 13
#data$fourteen<- 14
#data$fifteen<- 15
#data$sixteen<- 16

data$xvalue<- 0


#givine values to each nucecotide and if they have a cetain combo

for (i in 1:length(data$num)) {
    if (data$bigAAChange[i] == 0 && data$makesCpG[i] == 0 && data$WTnt[i] == "a") {
        data$xvalue[i] <- 1
    } else if (data$bigAAChange[i] == 1 && data$makesCpG[i] == 0 && data$WTnt[i] == "a") {
        data$xvalue[i] <- 2
    } else if (data$bigAAChange[i] == 0 && data$makesCpG[i] == 1 && data$WTnt[i] == "a") {
        data$xvalue[i] <- 3
    } else if (data$bigAAChange[i] == 1 && data$makesCpG[i] == 1 && data$WTnt[i] == "a") {
        data$xvalue[i] <- 4
    } else if (data$bigAAChange[i] == 0 && data$makesCpG[i] == 0 && data$WTnt[i] == "g") {
        data$xvalue[i] <- 5
    } else if (data$bigAAChange[i] == 1 && data$makesCpG[i] == 0 && data$WTnt[i] == "g") {
        data$xvalue[i] <- 6
    } else if (data$bigAAChange[i] == 0 && data$makesCpG[i] == 1 && data$WTnt[i] == "g") {
        data$xvalue[i] <- 7
    } else if (data$bigAAChange[i] == 1 && data$makesCpG[i] == 1 && data$WTnt[i] == "g") {
        data$xvalue[i] <- 8
    } else if (data$bigAAChange[i] == 0 && data$makesCpG[i] == 0 && data$WTnt[i] == "c") {
        data$xvalue[i] <- 9
    } else if (data$bigAAChange[i] == 1 && data$makesCpG[i] == 0 && data$WTnt[i] == "c") {
        data$xvalue[i] <- 10
    } else if (data$bigAAChange[i] == 0 && data$makesCpG[i] == 1 && data$WTnt[i] == "c") {
        data$xvalue[i] <- 11
    } else if (data$bigAAChange[i] == 1 && data$makesCpG[i] == 1 && data$WTnt[i] == "c") {
        data$xvalue[i] <- 12
    } else if (data$bigAAChange[i] == 0 && data$makesCpG[i] == 0 && data$WTnt[i] == "t") {
        data$xvalue[i] <- 13
    } else if (data$bigAAChange[i] == 1 && data$makesCpG[i] == 0 && data$WTnt[i] == "t") {
        data$xvalue[i] <- 14
    } else if (data$bigAAChange[i] == 0 && data$makesCpG[i] == 1 && data$WTnt[i] == "t") {
        data$xvalue[i] <- 15
    } else if (data$bigAAChange[i] == 1 && data$makesCpG[i] == 1 && data$WTnt[i] == "t") {
        data$xvalue[i] <- 16
    }
}

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
synYNa <- subset(synYNdata, WTnt=="a")
synYNc <- subset(synYNdata, WTnt=="c")
synYNg  <- subset(synYNdata, WTnt=="g")
synYNt  <- subset(synYNdata, WTnt=="t")
#syn subset for yes AA and yes CPG
synYYdata <- subset(syndata, combo=="yesAA yesCPG")
#syn subset for yes AA and yes CPG for a, c, g, t (none for HIV)
synYYa <- subset(synYYdata, WTnt=="a")
synYYc <- subset(synYYdata, WTnt=="c")
synYYg  <- subset(synYYdata, WTnt=="g")
synYYt  <- subset(synYYdata, WTnt=="t")


#graph 1 Synonymous Sites




#still working
#graph <- 
meansyn <- ((syndata$lowerConf+syndata$upperConf)/2)
ggplot(aes(factor(WTnt), MeanFreq, colour = combo), data = syndata)+
    #geom_jitter()+
    geom_point()+
    geom_errorbar(aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
    
#not yet done
graph + geom_boxplot() + 
    geom_dotplot(binaxis='y', 
                 stackdir='center', 
                 dotsize = .5, 
                 fill="blue") +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) 






#building the subset plot
ggplot(aes(factor(WTnt), MeanFreq), data = synNYa)+
    geom_jitter(col = "red") +
    geom_errorbar(aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(aes('a',median(c(median(lowerConf),median(upperConf)))))
# need to find lower and upper conf for each individual nuc (a,c,g,t)

ggplot(aes(factor(WTnt), MeanFreq), data = cpGdata)+
    geom_jitter(color ="blue")

ggplot()+
    #geom_jitter(data = synNYa, aes(factor(WTnt), MeanFreq),col = "red") +
    geom_errorbar(data = synNYa, aes(ymin = median(lowerConf), ymax = median(upperConf),position = 1, width = 0.2))+
    geom_point(aes('a',median(c(median(lowerConf),median(upperConf)))),data = synNYa)
?geom_errorbar


#created 10/11
ggplot(aes(factor(WTnt), MeanFreq), data = synNNdata)+
    geom_jitter(fill=5, col = "red") +
    geom_errorbar(data = alla, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =alla, aes('a',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = allc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =allc, aes('c',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = allg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =allg, aes('g',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = allt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =allt, aes('t',median(c(median(lowerConf),median(upperConf)))))


# combinging the differt colors
#In HIV data there is no C,G change for synNY values
ggplot(aes(factor(WTnt), MeanFreq), data = syndata)+
    #synNNdata
    geom_jitter(data = synNNdata, aes(factor(WTnt), MeanFreq),fill=5, col = "red", position=position_nudge(x=0.2), width=0.5) +
    geom_errorbar(data = synNNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNa, aes('a',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNc, aes('c',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNg, aes('g',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNt, aes('t',median(c(median(lowerConf),median(upperConf)))))+
    #synNYdata
    geom_jitter(data = synNYdata, aes(factor(WTnt), MeanFreq),fill=5, col = "blue", position=position_nudge(x=-0.2), width=0.5) +
    geom_errorbar(data = synNYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNYa, aes('a',median(c(median(lowerConf),median(upperConf)))))+
    #geom_errorbar(data = synNYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    #geom_point(data =synNYc, aes('c',median(c(median(lowerConf),median(upperConf)))))+
   # geom_errorbar(data = synNYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
   # geom_point(data =synNYg, aes('g',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNYt, aes('t',median(c(median(lowerConf),median(upperConf)))))

#In HIV data there is no C,G change for synNY values
ggplot(aes(factor(WTnt), MeanFreq), data = syndata)+
    #synNNdata
    geom_jitter(data = synNNdata, aes(factor(WTnt), MeanFreq),fill=5, col = "red", position=position_jitter(width = .05)) +
    geom_errorbar(data = synNNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNa, aes('a',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNc, aes('c',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNg, aes('g',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNt, aes('t',median(c(median(lowerConf),median(upperConf)))))+
    #synNYdata
    geom_jitter(data = synNYdata, aes(factor(WTnt), MeanFreq),fill=5, col = "blue", position=position_jitter(width = .05)) +
    geom_errorbar(data = synNYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNYa, aes('a',median(c(median(lowerConf),median(upperConf)))))+
    #geom_errorbar(data = synNYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    #geom_point(data =synNYc, aes('c',median(c(median(lowerConf),median(upperConf)))))+
    # geom_errorbar(data = synNYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    # geom_point(data =synNYg, aes('g',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNYt, aes('t',median(c(median(lowerConf),median(upperConf)))))



#new functions
ggplot(aes(factor(WTnt), MeanFreq), data = syndata)+
    #synNNdata
    scale_x_discrete() +
    geom_jitter(aes(colour = colsyn, x = factor(WTnt)),position = position_jitter(width = .05), alpha = 0.5) +
    #facet_wrap(~ factor(WTnt))+
    geom_errorbar(data = synNNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_boxplot(outlier.colour = NA, position = position_dodge(width=0.9))+
    geom_point(data =synNNa, aes('a',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNc, aes('c',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNg, aes('g',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNt, aes('t',median(c(median(lowerConf),median(upperConf)))))+
    #synNYdata
    geom_errorbar(data = synNYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNYa, aes('a',median(c(median(lowerConf),median(upperConf)))))+
    #geom_errorbar(data = synNYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    #geom_point(data =synNYc, aes('c',median(c(median(lowerConf),median(upperConf)))))+
    # geom_errorbar(data = synNYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    # geom_point(data =synNYg, aes('g',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNYt, aes('t',median(c(median(lowerConf),median(upperConf)))))


?facet_wrap()


#tring out numbers
ggplot(aes(factor(number), MeanFreq), data = syndata)+
    #synNNdata
    geom_jitter(data = synNNdata, aes(factor(number), MeanFreq),fill=5, col = "red") +
    geom_errorbar(data = synNNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNa, aes('1',median(c(median(lowerConf),median(upperConf)))))+
    geom_jitter(data = synNYdata, aes(factor(number), MeanFreq),fill=5, col = "blue")
    



    geom_jitter(col = "red") +
    geom_errorbar(aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(aes('a',median(c(median(lowerConf),median(upperConf)))))

    
    
    
    
#plotting with new xvalue 10/15
    ggplot(aes(factor(xvalue), MeanFreq), data = data)+
        #synNNdata 
        geom_jitter(data = synNNdata, aes(factor(xvalue), MeanFreq),fill=5, col = "red") +
        scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16))+
        geom_errorbar(data = synNNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
        geom_point(data =synNNa, aes('1',median(c(median(lowerConf),median(upperConf)))))+
        geom_errorbar(data = synNNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
        geom_point(data =synNNc, aes('9',median(c(median(lowerConf),median(upperConf)))))+
        geom_errorbar(data = synNNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
        geom_point(data =synNNg, aes('5',median(c(median(lowerConf),median(upperConf)))))+
        geom_errorbar(data = synNNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
        geom_point(data =synNNt, aes('13',median(c(median(lowerConf),median(upperConf)))))+
        #synNYdata
        geom_jitter(data = synNYdata, aes(factor(xvalue), MeanFreq),fill=5, col = "blue")+
        geom_errorbar(data = synNYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
        geom_point(data =synNYa, aes('3',median(c(median(lowerConf),median(upperConf)))))+
        #geom_errorbar(data = synNYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
        #geom_point(data =synNYc, aes('c',median(c(median(lowerConf),median(upperConf)))))+
        # geom_errorbar(data = synNYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
        # geom_point(data =synNYg, aes('g',median(c(median(lowerConf),median(upperConf)))))+
        geom_errorbar(data = synNYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
        geom_point(data =synNYt, aes('15',median(c(median(lowerConf),median(upperConf)))))
    
    
    
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
