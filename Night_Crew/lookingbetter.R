#Cleaning up for final project

setwd("~/Desktop/Git")

# freq calc v2 GL

library(seqinr)
 
bk <- read.fasta("bk.txt")

# reference sequence
ref <- bk[[1]]

# dataframe columns
num <- c(1:1089)
wtnt <- c()
freq <- c()

# for freq calculation later
absfreq <- c(rep(0,1089))
totalcount <- c(rep(0,1089))

# average WT calculation
# counts number of each nucleotide in each position
acount <- c(rep(0,1089))
gcount <- c(rep(0,1089))
ccount <- c(rep(0,1089))
tcount <- c(rep(0,1089))
nuc <- c()

# same as line 20 comment
for (i in 1:length(bk)) {
    sequence <- bk[[i]]
    for (j in 1:length(sequence)) {
        if (sequence[j] == 'a') {acount[j] = acount[j] + 1}
        if (sequence[j] == 'g') {gcount[j] = gcount[j] + 1}
        if (sequence[j] == 'c') {ccount[j] = ccount[j] + 1}
        if (sequence[j] == 't') {tcount[j] = tcount[j] + 1}
    }
}

# assigns wtnt based on most frequent nucleotide across all sequences per position
for (j in 1:length(sequence)) {
    nuc[j] <- max(c(acount[j],gcount[j],ccount[j],tcount[j]))
    if (max(nuc[j]) == acount[j]) {wtnt[j] <- 'a'}
    if (max(nuc[j]) == gcount[j]) {wtnt[j] <- 'g'}
    if (max(nuc[j]) == ccount[j]) {wtnt[j] <- 'c'}
    if (max(nuc[j]) == tcount[j]) {wtnt[j] <- 't'}
    nuc <- c()
}

# gives absolute totals to be used for frequency calculation
for (i in 1:length(bk)) {
    sequence <- bk[[i]]
    for (j in 1:length(sequence)){
        if (wtnt[j] == 'a') {
            if (sequence[j] == 'g') {
                absfreq[j] <- absfreq[j] + 1
                totalcount[j] <- totalcount[j] + 1
            }
            else if (sequence[j] == 'a') {
                totalcount[j] <- totalcount[j] + 1
            }
        }
        if (wtnt[j] == 'g') {
            if (sequence[j] == 'a') {
                absfreq[j] <- absfreq[j] + 1
                totalcount[j] <- totalcount[j] + 1
            }
            else if (sequence[j] == 'g') {
                totalcount[j] <- totalcount[j] + 1
            }
        }
        if (wtnt[j] == 'c') {
            if (sequence[j] == 't') {
                absfreq[j] <- absfreq[j] + 1
                totalcount[j] <- totalcount[j] + 1
            }
            else if (sequence[j] == 'c') {
                totalcount[j] <- totalcount[j] + 1
            }
        }
        if (wtnt[j] == 't') {
            if (sequence[j] == 'c') {
                absfreq[j] <- absfreq[j] + 1
                totalcount[j] <- totalcount[j] + 1
            }
            else if (sequence[j] == 't') {
                totalcount[j] <- totalcount[j] + 1
            }
        }
    }
}

# calculates frequency as percentage
for (i in 1:length(absfreq)) {
    freq[i] <- absfreq[i] / totalcount[i]
}

# creates dataframe containing all data
bk_data <- data.frame(num,wtnt,freq)

# outputs data to file
write.csv(bk_data,"bk_data.csv")

#Start of the graph


#install.packages('ggplot2')
library(ggplot2)

#reading in data
data<-read.csv('OverviewSelCoeff_BachelerFilter.csv') 

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


#Subsetting for the data that we really want
datatww <- subset(data, TypeOfSite=="syn" | TypeOfSite=="nonsyn")


#givine values to each nucecotide and if they have a cetain combo
#first need to introduce new columes for collecting the xvaue and the color
datatww$xvalue<- 0
datatww$color <- 0

for (i in 1:length(datatww$num)) {
    if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "a") {
        datatww$xvalue[i] <- 1
        datatww$color[i] <- "blue"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "a") {
        datatww$xvalue[i] <- 2
        datatww$color[i] <- "red"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "a") {
        datatww$xvalue[i] <- 3
        datatww$color[i] <- "green"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "a") {
        datatww$xvalue[i] <- 4
        datatww$color[i] <- "purple"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "g") {
        datatww$xvalue[i] <- 5
        datatww$color[i] <- "blue"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "g") {
        datatww$xvalue[i] <- 6
        datatww$color[i] <- "red"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "g") {
        datatww$xvalue[i] <- 7
        datatww$color[i] <- "green"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "g") {
        datatww$xvalue[i] <- 8
        datatww$color[i] <- "purple"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "c") {
        datatww$xvalue[i] <- 9
        datatww$color[i] <- "blue"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "c") {
        datatww$xvalue[i] <- 10
        datatww$color[i] <- "red"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "c") {
        datatww$xvalue[i] <- 11
        datatww$color[i] <- "green"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "c") {
        datatww$xvalue[i] <- 12
        datatww$color[i] <- "purple"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "t") {
        datatww$xvalue[i] <- 13
        datatww$color[i] <- "blue"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "t") {
        datatww$xvalue[i] <- 14
        datatww$color[i] <- "red"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "t") {
        datatww$xvalue[i] <- 15
        datatww$color[i] <- "green"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "t") {
        datatww$xvalue[i] <- 16
        datatww$color[i] <- "purple"
    }
}


#subsetiing all the data by amino acid, drastic AA, and CpG forming 

#syn subset stuff
syndata <- subset(datatww, TypeOfSite=="syn")

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

#nonsyn sub set stuff
nonsyndata <- subset(datatww, TypeOfSite=="nonsyn")

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
nonsynYNa <- subset(nonYNdata, WTnt=="a")
nonsynYNc <- subset(nonYNdata, WTnt=="c")
nonsynYNg  <- subset(nonYNdata, WTnt=="g")
nonsynYNt  <- subset(nonYNdata, WTnt=="t")
#syn subset for yes AA and yes CPG
nonYYdata <- subset(nonsyndata, combo=="yesAA yesCPG")
#syn subset for yes AA and yes CPG for a, c, g, t (for HIV a, t)
nonsynYYa <- subset(nonYYdata, WTnt=="a")
nonsynYYc <- subset(nonYYdata, WTnt=="c")
nonsynYYg  <- subset(nonYYdata, WTnt=="g")
nonsynYYt  <- subset(nonYYdata, WTnt=="t")


#datatww$TypeOfSite <- revalue(x = datatww$TypeOfSite , c("syn" = "Synonymous Sites", "nonsyn" = "Non-synonymous Sites"))
#levels(datatww$TypeOfSite) <- gsub("syn", "Synonymous Sites", levels(datatww$TypeOfSite))
#levels(datatww$TypeOfSite) <- gsub("nonsyn", "Nonsynonymous Sites", levels(datatww$TypeOfSite))
#still need to log
#graph
ggplot(aes(factor(xvalue), MeanFreq), data = datatww)+
    #synNNdata
    scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),breaks=c("2","6","10", "14"), labels=c("a", "g", "c","t"))+
    geom_jitter(data= syndata,aes(colour = syndata$color, x = factor(xvalue)),position = position_jitter(width = .2), alpha = 0.5) +
    facet_wrap(~ TypeOfSite)+
    geom_errorbar(data = synNNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNa, aes('1',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNc, aes('9',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNg, aes('5',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNNt, aes('13',median(c(median(lowerConf),median(upperConf)))))+
    #synNYdata
    geom_errorbar(data = synNYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNYa, aes('3',median(c(median(lowerConf),median(upperConf)))))+
    #geom_errorbar(data = synNYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    #geom_point(data =synNYc, aes('10',median(c(median(lowerConf),median(upperConf)))))+
    # geom_errorbar(data = synNYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    # geom_point(data =synNYg, aes('6',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = synNYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =synNYt, aes('15',median(c(median(lowerConf),median(upperConf)))))+
    
    
    geom_jitter(data= nonsyndata,aes(colour = nonsyndata$color, x = factor(xvalue)),position = position_jitter(width = .2), alpha = 0.5) +
    geom_errorbar(data = nonsynNNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynNNa, aes('1',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = nonsynNNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynNNc, aes('9',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = nonsynNNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynNNg, aes('5',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = nonsynNNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynNNt, aes('13',median(c(median(lowerConf),median(upperConf)))))+
    #
    geom_errorbar(data = nonsynNYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynNYa, aes('3',median(c(median(lowerConf),median(upperConf)))))+
    #geom_errorbar(data = nonsynNYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    #geom_point(data =nonsynNYc, aes('10',median(c(median(lowerConf),median(upperConf)))))+
    #geom_errorbar(data = nonsynNYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    #geom_point(data =nonsynNYg, aes('6',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = nonsynNYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynNYt, aes('15',median(c(median(lowerConf),median(upperConf)))))+
    #
    geom_errorbar(data = nonsynYNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynYNa, aes('2',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = nonsynYNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynYNc, aes('10',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = nonsynYNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynYNg, aes('6',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = nonsynYNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynYNt, aes('14',median(c(median(lowerConf),median(upperConf)))))+
    #
    geom_errorbar(data = nonsynYYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynYYa, aes('4',median(c(median(lowerConf),median(upperConf)))))+
    #geom_errorbar(data = nonsynYYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    #geom_point(data =nonsynYYc, aes('10',median(c(median(lowerConf),median(upperConf)))))+
    #geom_errorbar(data = nonsynYYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    #geom_point(data =nonsynYYg, aes('6',median(c(median(lowerConf),median(upperConf)))))+
    geom_errorbar(data = nonsynYYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))+
    geom_point(data =nonsynYYt, aes('16',median(c(median(lowerConf),median(upperConf)))))+
    
    scale_color_manual(labels = c("No drastic AA change (non-Cpg-forming)","No drastic AA change (Cpg-forming)","Drastic AA change (non-Cpg-forming)","Drastic AA change (Cpg-forming)"), values = c("blue", "red","green", "purple")) +
    labs(x="Mutation Type", y="Mutation Frquency",  
         col=" ")
    