#Cleaning up for final project

setwd("~/Desktop/Git")

#install.packages('ggplot2')
#install.packages('scales')
#scales is used to log the y axis in the graoh and ggplot2 is used to create the graph
library(scales)
library(ggplot2)

#reading in data
data<-read.csv('OverviewSelCoeff_BachelerFilter.csv') 

#Start of the graph function
nightcrewgraph = function(data){
#building the combo lines to help sort (used later when info is subsetted) 
data$combo<- (data$bigAAChange*3) + (data$makesCpG*2)
data$combo<- as.factor(data$combo)
levels(data$combo) <- gsub("0", "noAA noCPG", levels(data$combo))
levels(data$combo) <- gsub("2", "noAA yesCPG", levels(data$combo))
levels(data$combo) <- gsub("3", "yesAA noCPG", levels(data$combo))
levels(data$combo) <- gsub("5", "yesAA yesCPG", levels(data$combo))

#Subsetting for the data that we really want
datatww <- subset(data, TypeOfSite=="syn" | TypeOfSite=="nonsyn")
datatww$TypeOfSite <- ifelse(datatww$TypeOfSite == "syn", "Synonymous Sites", "Non-synonymous Sites")
#idea for log
#datatww$lowerConf= log(datatww$lowerConf)
#datatww$upperConf= log(datatww$upperConf)
#datatww$MeanFreq = log(datatww$MeanFreq) 

#givine values to each nucecotide and if they have a cetain combo
#first need to introduce new columes for collecting the xvaue and the color
datatww$xvalue<- 0
datatww$color <- 0

for (i in 1:length(datatww$num)) {
    if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "a") {
        datatww$xvalue[i] <- 1
        datatww$color[i] <- "blue"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "a") {
        datatww$xvalue[i] <- 2
        datatww$color[i] <- "red"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "a") {
        datatww$xvalue[i] <- 3
        datatww$color[i] <- "green"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "a") {
        datatww$xvalue[i] <- 4
        datatww$color[i] <- "purple"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "g") {
        datatww$xvalue[i] <- 5
        datatww$color[i] <- "blue"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "g") {
        datatww$xvalue[i] <- 6
        datatww$color[i] <- "red"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "g") {
        datatww$xvalue[i] <- 7
        datatww$color[i] <- "green"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "g") {
        datatww$xvalue[i] <- 8
        datatww$color[i] <- "purple"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "c") {
        datatww$xvalue[i] <- 9
        datatww$color[i] <- "blue"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "c") {
        datatww$xvalue[i] <- 10
        datatww$color[i] <- "red"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "c") {
        datatww$xvalue[i] <- 11
        datatww$color[i] <- "green"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "c") {
        datatww$xvalue[i] <- 12
        datatww$color[i] <- "purple"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "t") {
        datatww$xvalue[i] <- 13
        datatww$color[i] <- "blue"
    } else if (datatww$bigAAChange[i] == 0 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "t") {
        datatww$xvalue[i] <- 14
        datatww$color[i] <- "red"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 0 && datatww$WTnt[i] == "t") {
        datatww$xvalue[i] <- 15
        datatww$color[i] <- "green"
    } else if (datatww$bigAAChange[i] == 1 && datatww$makesCpG[i] == 1 && datatww$WTnt[i] == "t") {
        datatww$xvalue[i] <- 16
        datatww$color[i] <- "purple"
    }
}




#subsetiing all the data by amino acid, drastic AA, and CpG forming 

#syn subset stuff
syndata <- subset(datatww, TypeOfSite=="Synonymous Sites")

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
nonsyndata <- subset(datatww, TypeOfSite=="Non-synonymous Sites")

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


#graph
nightmare <-ggplot(aes(factor(xvalue), MeanFreq), data = datatww)+
    #synNNdata
    scale_y_log10() +
    scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),breaks=c("2","6","10", "14"), labels=c("a -> g", "g -> a", "c -> t","t -> c"))+
    geom_jitter(data= syndata,aes(colour = syndata$color, x = factor(xvalue)),position = position_jitter(width = .2), alpha = 0.5) +
    facet_wrap(~ TypeOfSite)+
    #nonsyn data
    geom_jitter(data= nonsyndata,aes(colour = nonsyndata$color, x = factor(xvalue)),position = position_jitter(width = .2), alpha = 0.5) +
    geom_vline(xintercept = 4.5, linetype="dotted", color = "dark gray", size=1.5) +
    geom_vline(xintercept = 8.1, linetype="dotted", color = "dark gray", size=1.5) +
    geom_vline(xintercept = 12, linetype="dotted", color = "dark gray", size=1.5)+
    scale_color_manual(labels = c("No drastic AA change (non-Cpg-forming)","Drastic AA change (non-Cpg-forming)","Drastic AA change (Cpg-forming)", "No drastic AA change (Cpg-forming)"), values = c("green", "yellow","red", "blue")) +
    labs(x="Mutation Type", y="Mutation Frquency",col=" ")

#graphing other info basied on conditions
nightmare
    if (nrow(synNNa)!=0) {
        nightmare<- nightmare + geom_errorbar(data = synNNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare + geom_point(data =synNNa, aes('1',median(c(median(lowerConf),median(upperConf)))))
    }
    if (nrow(synNNc)!=0) {
        nightmare<- nightmare + geom_errorbar(data = synNNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare + geom_point(data =synNNc, aes('9',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synNNg)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synNNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synNNg, aes('5',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synNNt)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synNNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synNNt, aes('13',median(c(median(lowerConf),median(upperConf)))))
    } 
    #synNYdata
    if (nrow(synNYa)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synNYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synNYa, aes('2',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synNYc)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synNYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synNYc, aes('10',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synNYg)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synNYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synNYg, aes('6',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synNYt)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synNYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synNYt, aes('14',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synYNa)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synYNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synYNa, aes('3',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synYNc)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synYNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synYNc, aes('11',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synYNg)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synYNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synYNg, aes('7',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synYNt)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synYNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synYNt, aes('15',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synYYa)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synYYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synYYa, aes('4',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synYYc)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synYYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synYYc, aes('12',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synYYg)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synYYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synYYg, aes('8',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(synYYt)!=0) {
        nightmare<- nightmare +geom_errorbar(data = synYYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =synYYt, aes('16',median(c(median(lowerConf),median(upperConf)))))
    } 
    
    #nonsyn data
    
    if (nrow(nonsynNNa)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynNNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynNNa, aes('1',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynNNc)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynNNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynNNc, aes('9',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynNNg)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynNNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynNNg, aes('5',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynNNt)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynNNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynNNt, aes('13',median(c(median(lowerConf),median(upperConf)))))
    } 
    
    if (nrow(nonsynNYa)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynNYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare + geom_point(data =nonsynNYa, aes('2',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynNYc)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynNYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynNYc, aes('10',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynNYg)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynNYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynNYg, aes('6',median(c(median(lowerConf),median(upperConf)))))
    }  
    if (nrow(nonsynNYt)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynNYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynNYt, aes('14',median(c(median(lowerConf),median(upperConf)))))
    } 
    #
    if (nrow(nonsynYNa)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynYNa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynYNa, aes('3',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynYNc)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynYNc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynYNc, aes('11',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynYNg)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynYNg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynYNg, aes('7',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynYNt)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynYNt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynYNt, aes('15',median(c(median(lowerConf),median(upperConf)))))
    } 
    #
    if (nrow(nonsynYYa)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynYYa, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynYYa, aes('4',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynYYc)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynYYc, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynYYc, aes('12',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynYYg)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynYYg, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynYYg, aes('8',median(c(median(lowerConf),median(upperConf)))))
    } 
    if (nrow(nonsynYYt)!=0) {
        nightmare<- nightmare +geom_errorbar(data = nonsynYYt, aes(ymin = median(lowerConf), ymax = median(upperConf), width = 0.2))
        nightmare<- nightmare +geom_point(data =nonsynYYt, aes('16',median(c(median(lowerConf),median(upperConf)))))
    } 
nightmare <- nightmare +  theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank()) 


print(nightmare)

#return(datatww)
#ggsave(filename="nightmare.pdf", plot=nightmare)
#ggsave(filename="nightmare.pdf", plot = last_plot(), device = NULL, path = NULL,
       #scale = 1, width = NA, height = NA, units = c("in", "cm", "mm"),
       #dpi = 300, limitsize = TRUE, ...)
}

nightcrewgraph(data)
