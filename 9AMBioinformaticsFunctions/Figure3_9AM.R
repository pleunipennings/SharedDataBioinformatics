#Night_Crew Team
#Members Victoria, Sam, Gordon, Mordy 
#setwd("~/Desktop/Git")

#install.packages('ggplot2')
#install.packages('scales')
#scales is used to log the y axis in the graoh and ggplot2 is used to create the graph
library(scales)
library(ggplot2)

#this function graphs transition mutations separated by Synonymous Sites and Nonsynonymous sites
# this graph show the frequencies of the changes and also seperates combinations of drastic AA and CpG sites by color
#Start of the graph function, here data is your dataset. Make sure you have these columns and the column names are the same.bigAAChange, makesCpG, TypeOfSite, num, MeanFreqor or freq, wtnt or WTnt 
#In column TypeOfSite make sure that syn, nonsyn are present. No - or full names
comparing_mutation_graph = function(data){
    #building the combo lines to help sort (used later when info is subsetted) 
    data$combo<- (data$bigAAChange*3) + (data$makesCpG*2)
    data$combo<- as.factor(data$combo)
    levels(data$combo) <- gsub("0", "noAA noCPG", levels(data$combo))
    levels(data$combo) <- gsub("2", "noAA yesCPG", levels(data$combo))
    levels(data$combo) <- gsub("3", "yesAA noCPG", levels(data$combo))
    levels(data$combo) <- gsub("5", "yesAA yesCPG", levels(data$combo))
    #This will change your colomn name to what we orginally used to make our graph
    colnames(data)[which(names(data) == "wtnt")] <- "WTnt"
    colnames(data)[which(names(data) == "freq")] <- "MeanFreq"
    #Subsetting for the data that we really want just Synonymous and nonynonymous
    datatww <- subset(data, TypeOfSite=="syn" | TypeOfSite=="nonsyn")
    #if your dataframe has "-" try using the following line of code
    #datatww<-datatww[!(which(datatww$WTnt=="-")),]
    datatww$TypeOfSite <- ifelse(datatww$TypeOfSite == "syn", "Synonymous Sites", "Non-synonymous Sites")
    
    #gives values to each nucecotide and if they have a cetain combo
    #first need to introduce new columes for collecting the xvalue and the color
    datatww$xvalue<- 0 #giving each nucecotide a value depending on whether or not there is an AA change a CpG site and the nucecotide
    datatww$color <- 0 # four colors for the four combos: Blue = no bigAAChange/no makesCpG, Red = no bigAAChange/yes makesCpG, Green = yes bigAAChange/no makesCpG, Purple = yes bigAAChange/yes makesCpG 
    #for loop adds color and a xvalue to the dataset, for all nucecotides
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
    #subsetiing all the data by amino acid, drastic AA, and CpG forming: Helps with analyst 
    #syn subset 
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
    
    #nonsyn sub set 
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
    
    
    #the actual graph is here.
    graph <-ggplot(aes(factor(xvalue), MeanFreq), data = datatww)+
        #log scale to make the data eaisier to see
        scale_y_log10() +
        #scale_x_discrete makes a point at each number, breaks lets us seperate into sections for each nucecotide
        scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),breaks=c("2","6","10", "14"), labels=c("a -> g", "g -> a", "c -> t","t -> c"))+
        #syn data, jitter seperates the points
        geom_jitter(data= syndata,aes(colour = syndata$color, x = factor(xvalue)),position = position_jitter(width = .2), alpha = 0.5) +
        #facet_wrap splits graph between syn and nonsyn
        facet_wrap(~ TypeOfSite)+
        #nonsyn data, jitter seperates the points
        geom_jitter(data= nonsyndata,aes(colour = nonsyndata$color, x = factor(xvalue)),position = position_jitter(width = .2), alpha = 0.5) +
        #draws lines between the nucelotides 
        geom_vline(xintercept = 4.5, linetype="solid", color = "black", size=0.9) +
        geom_vline(xintercept = 8.1, linetype="solid", color = "black", size=0.9) +
        geom_vline(xintercept = 12, linetype="solid", color = "black", size=0.9)+
        #give points new colors and lables the colors
        scale_color_manual(labels = c("No drastic AA change (non-Cpg-forming)","Drastic AA change (non-Cpg-forming)","Drastic AA change (Cpg-forming)", "No drastic AA change (Cpg-forming)"), values = c("firebrick", "darkolivegreen","goldenrod3", "royalblue3")) +
        #labels X and Y axis
        labs(x="Mutation Type", y="Mutation Frquency",col=" ")
    
    #graphing other info basied on conditions, used as a check to see if data is being read corectly
    #if statments checks to see if the dataset exist. If things that should exist show up on the graph you would be able to determine what subset of data needs to be checked
    if (nrow(synNNa)!=0) {
        synNNaconf <- t.test(synNNa$MeanFreq)$conf.int
        graph<- graph + geom_errorbar(data = synNNa, aes(ymin = synNNaconf[[1]], ymax = synNNaconf[[2]], width = 0.2))
        graph<- graph + geom_point(data =synNNa, aes('1',mean(synNNa$MeanFreq)))
    }
    if (nrow(synNNc)!=0) {
        synNNcconf <- t.test(synNNc$MeanFreq)$conf.int
        graph<- graph + geom_errorbar(data = synNNc, aes(ymin = synNNcconf[[1]], ymax = synNNcconf[[2]], width = 0.2))
        graph<- graph + geom_point(data =synNNc, aes('9',mean(synNNc$MeanFreq)))
    } 
    if (nrow(synNNg)!=0) {
        synNNgconf <- t.test(synNNg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNNg, aes(ymin = synNNgconf[[1]], ymax = synNNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNNg, aes('5',mean(synNNg$MeanFreq)))
    } 
    if (nrow(synNNt)!=0) {
        synNNtconf <- t.test(synNNt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNNt, aes(ymin = synNNtconf[[1]], ymax = synNNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNNt, aes('13',mean(synNNt$MeanFreq)))
    } 
    #synNYdata
    if (nrow(synNYa)!=0) {
        synNYaconf <- t.test(synNYa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNYa, aes(ymin = synNYaconf[[1]], ymax = synNYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYa, aes('2',mean(synNYa$MeanFreq)))
    } 
    if (nrow(synNYc)!=0) {
        synNYcconf <- t.test(synNYc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNYc, aes(ymin = synNYcconf[[1]], ymax = synNYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYc, aes('10',mean(synNYc$MeanFreq)))
    } 
    if (nrow(synNYg)!=0) {
        synNYgconf <- t.test(synNYg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNYg, aes(ymin = synNYgconf[[2]], ymax = synNYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYg, aes('6',mean(synNYg$MeanFreq)))
    } 
    if (nrow(synNYt)!=0) {
        synNYtconf <- t.test(synNYt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synNYt, aes(ymin = synNYtconf[[1]], ymax = synNYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synNYt, aes('14',mean(synNYt$MeanFreq)))
    } 
    if (nrow(synYNa)!=0) {
        synYNaconf <- t.test(synYNa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYNa, aes(ymin = synYNaconf[[1]], ymax = synYNaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNa, aes('3',mean(synYNa$MeanFreq)))
    } 
    if (nrow(synYNc)!=0) {
        synYNcconf <- t.test(synYNc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYNc, aes(ymin = synYNcconf[[1]], ymax = synYNcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNc, aes('11',mean(synYNc$MeanFreq)))
    } 
    if (nrow(synYNg)!=0) {
        synYNgconf <- t.test(synYNg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYNg, aes(ymin = synYNgconf[[1]], ymax = synYNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNg, aes('7',mean(synYNg$MeanFreq)))
    } 
    if (nrow(synYNt)!=0) {
        synYNtconf <- t.test(synYNt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYNt, aes(ymin = synYNtconf[[1]], ymax = synYNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYNt, aes('15',mean(synYNt$MeanFreq)))
    } 
    if (nrow(synYYa)!=0) {
        synYYaconf <- t.test(synYYa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYYa, aes(ymin = synYYaconf[[1]], ymax = synYYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYa, aes('4',mean(synYYa$MeanFreq)))
    } 
    if (nrow(synYYc)!=0) {
        synYYcconf <- t.test(synYYc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYYc, aes(ymin = synYYcconf[[1]], ymax = synYYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYc, aes('12', mean(synYYc$MeanFreq)))
    } 
    if (nrow(synYYg)!=0) {
        synYYgconf <- t.test(synYYg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYYg, aes(ymin = synYYgconf[[1]], ymax = synYYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYg, aes('8',mean(synYYg$MeanFreq)))
    } 
    if (nrow(synYYt)!=0) {
        synYYtconf <- t.test(synYYt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = synYYt, aes(ymin = synYYtconf[[1]], ymax = synYYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =synYYt, aes('16',mean(synYYt$MeanFreq)))
    }
    #nonsyn data
    if (nrow(nonsynNNa)!=0) {
        nonsynNNaconf <- t.test(nonsynNNa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNa, aes(ymin = nonsynNNaconf[[1]], ymax = nonsynNNaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNa, aes('1',mean(nonsynNNa$MeanFreq)))
    } 
    if (nrow(nonsynNNc)!=0) {
        nonsynNNcconf <- t.test(nonsynNNc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNc, aes(ymin = nonsynNNcconf[[1]], ymax = nonsynNNcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNc, aes('9',mean(nonsynNNc$MeanFreq)))
    } 
    if (nrow(nonsynNNg)!=0) {
        nonsynNNgconf <- t.test(nonsynNNg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNg, aes(ymin = nonsynNNgconf[[1]], ymax = nonsynNNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNg, aes('5',mean(nonsynNNg$MeanFreq)))
    } 
    if (nrow(nonsynNNt)!=0) {
        nonsynNNtconf <- t.test(nonsynNNt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNNt, aes(ymin = nonsynNNtconf[[1]], ymax = nonsynNNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNNt, aes('13',mean(nonsynNNt$MeanFreq)))
    } 
    if (nrow(nonsynNYa)!=0) {
        nonsynNYaconf <- t.test(nonsynNYa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYa, aes(ymin = nonsynNYaconf[[1]], ymax = nonsynNYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYa, aes('2',mean(nonsynNYa$MeanFreq)))
    } 
    if (nrow(nonsynNYc)!=0) {
        nonsynNYcconf <- t.test(nonsynNYc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYc, aes(ymin = nonsynNYcconf[[1]], ymax = nonsynNYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYc, aes('10',mean(nonsynNYc$MeanFreq)))
    } 
    if (nrow(nonsynNYg)!=0) {
        nonsynNYgconf <- t.test(nonsynNYg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYg, aes(ymin = nonsynNYgconf[[1]], ymax = nonsynNYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYg, aes('6',mean(nonsynNYg$MeanFreq)))
    }  
    if (nrow(nonsynNYt)!=0) {
        nonsynNYtconf <- t.test(nonsynNYt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynNYt, aes(ymin = nonsynNYtconf[[1]], ymax = nonsynNYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynNYt, aes('14',mean(nonsynNYt$MeanFreq)))
    } 
    if (nrow(nonsynYNa)!=0) {
        nonsynYNaconf <- t.test(nonsynYNa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNa, aes(ymin = nonsynYNaconf[[1]], ymax = nonsynYNaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNa, aes('3',mean(nonsynYNa$MeanFreq)))
    } 
    if (nrow(nonsynYNc)!=0) {
        nonsynYNcconf <- t.test(nonsynYNc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNc, aes(ymin = nonsynYNcconf[[1]], ymax = nonsynYNcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNc, aes('11',mean(nonsynYNc$MeanFreq)))
    } 
    if (nrow(nonsynYNg)!=0) {
        nonsynYNgconf <- t.test(nonsynYNg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNg, aes(ymin = nonsynYNgconf[[1]], ymax = nonsynYNgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNg, aes('7',mean(nonsynYNg$MeanFreq)))
    } 
    if (nrow(nonsynYNt)!=0) {
        nonsynYNtconf <- t.test(nonsynYNt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYNt, aes(ymin = nonsynYNgconf[[1]], ymax = nonsynYNtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYNt, aes('15',mean(nonsynYNt$MeanFreq)))
    } 
    if (nrow(nonsynYYa)!=0) {
        nonsynYYaconf <- t.test(nonsynYYa$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYa, aes(ymin = nonsynYYaconf[[1]], ymax = nonsynYYaconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYa, aes('4',mean(nonsynYYa$MeanFreq)))
    } 
    if (nrow(nonsynYYc)!=0) {
        nonsynYYcconf <- t.test(nonsynYYc$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYc, aes(ymin = nonsynYYcconf[[1]], ymax = nonsynYYcconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYc, aes('12',mean(nonsynYYc$MeanFreq)))
    } 
    if (nrow(nonsynYYg)!=0) {
        nonsynYYgconf <- t.test(nonsynYYg$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYg, aes(ymin = nonsynYYgconf[[1]], ymax = nonsynYYgconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYg, aes('8',mean(nonsynYYg$MeanFreq)))
    } 
    if (nrow(nonsynYYt)!=0) {
        nonsynYYtconf <- t.test(nonsynYYt$MeanFreq)$conf.int
        graph<- graph +geom_errorbar(data = nonsynYYt, aes(ymin = nonsynYYtconf[[1]], ymax = nonsynYYtconf[[2]], width = 0.2))
        graph<- graph +geom_point(data =nonsynYYt, aes('16',mean(nonsynYYt$MeanFreq)))
    } 
    #this get rid of some of the back ground lines
    graph <- graph +  theme_bw() +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank(),
              #panel.border = element_blank(),
              panel.background = element_blank()) 
    
    
    print(graph)
}

comparing_mutation_graph(bk_data)
#data is your dataframe



