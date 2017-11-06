Fig3<-function(dfx){
  #begin plot function for figure 3; ends on line 651
 # df<-dffin3
  library(ggplot2)
  library(plyr)
  library(grid)
  library(scales)
  library(gridExtra)
  CSV<-dfx
  
  
  #This portion establishes plotting order => ggplot will plot things in aplhabetic order. The selections of data from the incoming data.frame will be tagged a thru r.
  
  ##SYN SITES (LEFT GRAPH)
  #all green points for the left synonomous site grapha
  
  ##SYN SITES (LEFT GRAPH)
  #all green points for the left synonomous site graphs
  a <- frequenciesOfSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "syn"  )) & (CSV$WTnt == "a" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
  summary(a)
  c <- frequenciesOfSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "syn"  )) & (CSV$WTnt == "t" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
  e <- frequenciesOfSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "syn"  )) & (CSV$WTnt == "c" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
  f <- frequenciesOfSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "syn"  )) & (CSV$WTnt == "g" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
  
  #blue dots of left graph that show a CpG-forming mutation
  b <- frequenciesOfSynAmutsCP <- CSV[which(((CSV$TypeOfSite == "syn"  )) & (CSV$WTnt == "a") & (CSV$bigAAChange == "0") & (CSV$makesCpG == "1")),"MeanFreq"]
  d <- frequenciesOfSynTmutsCP <- CSV[which(((CSV$TypeOfSite == "syn"  )) & (CSV$WTnt == "t") & (CSV$bigAAChange == "0") & (CSV$makesCpG == "1")),"MeanFreq"]
  
  #NON SYNONYMOUS SITES (RIGHT GRAPH)
  #all green points for the right nonsynonomous site graph
  g <- frequenciesOfNONSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "a" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
  k <- frequenciesOfNONSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "t" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
  o <- frequenciesOfNONSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "c" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
  q <- frequenciesOfNONSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "g" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
  
  #blue dots of right graph that show a CpG-forming mutation
  h <- frequenciesOfNONSynAmutsCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "a") & (CSV$bigAAChange == "0") & (CSV$makesCpG == "1")),"MeanFreq"]
  l <- frequenciesOfNONSynTmutsCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "t") & (CSV$bigAAChange == "0") & (CSV$makesCpG == "1")),"MeanFreq"]
  
  #all yellow points for the right nonsyn site graph
  i <- frequenciesOfNONSynAmutsDRASTIC <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "a" & (CSV$bigAAChange == "1") &(CSV$makesCpG == "0"))),"MeanFreq"]
  m <- frequenciesOfNONSynTmutsDRASTIC <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "t" & (CSV$bigAAChange == "1") &(CSV$makesCpG == "0"))),"MeanFreq"]
  p <- frequenciesOfNONSynGmutsDRASTIC <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "c" & (CSV$bigAAChange == "1") &(CSV$makesCpG == "0"))),"MeanFreq"]
  r <- frequenciesOfNONSynCmutsDRASTIC <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "g" & (CSV$bigAAChange == "1") &(CSV$makesCpG == "0"))),"MeanFreq"]
  
  
  #red points on the right
  j <- frequenciesOfNONSynAmutsCPDRASTIC <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "a" & (CSV$bigAAChange == "1") &(CSV$makesCpG == "1"))),"MeanFreq"]
  n <- frequenciesOfNONSynTmutsCPDRASTIC <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "t" & (CSV$bigAAChange == "1") &(CSV$makesCpG == "1"))),"MeanFreq"]
  
  #Since the results of the data selections are all vectors of differing lengths, they are collected into a list with tags according to plotting order
  mylist <- list (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) 
  k
  mylist
  namvec<-c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r") 
  names(mylist)<-namvec
  (mylist)
  #The list has some NULL elements; these must be replaced with NA since the rest of the code requires values in all lettered locations
  mylistvec<-c()
  for (i in 1:length(mylist)){
    if (length(mylist[[i]])==0){mylist[[i]] = NA}
    mylistvec<-c(mylistvec,mylist[i])
  }
  
  #ggplot requires that data be in a data.frame and to make most effective use of the plotting functions, data will be arranged in to a column of values 
  #and a column of plotting order identifiers (refvec)
  dfmylist<-data.frame(mylistvec[1])
  dfmylist$refvec<-names(mylistvec[1])
  colnames(dfmylist)<-c("values","refvec")
  for (i in 2:length(mylist)){
    dfwkg<-data.frame((mylistvec[i]))
    dfwkg$refvec<-names(mylistvec[i])
    colnames(dfwkg)<-c("values","refvec")
    dfmylist<-rbind(dfmylist,dfwkg)
  }
  #View (dfmylist)
  
  
  
  #Now add a column for AA_Category and a column for colors so tha ggplot can pull colors and legend from data.frame 
  dfmylist[which(dfmylist$refvec == "a"),"AA_category"]<-"Non-drastic, non-CpG" 
  dfmylist[which(dfmylist$refvec == "b"),"AA_category"]<-"Non-drastic, CpG" 
  dfmylist[which(dfmylist$refvec == "a"),"AA_category"]<-"Non-drastic, non-CpG" 
  dfmylist[which(dfmylist$refvec == "b"),"AA_category"]<-"Non-drastic, CpG" 
  dfmylist[which(dfmylist$refvec == "c"),"AA_category"]<-"Non-drastic, non-CpG" 
  dfmylist[which(dfmylist$refvec == "d"),"AA_category"]<-"Non-drastic, CpG" 
  dfmylist[which(dfmylist$refvec == "e"),"AA_category"]<-"Non-drastic, non-CpG" 
  dfmylist[which(dfmylist$refvec == "f"),"AA_category"]<-"Non-drastic, non-CpG" 
  
  dfmylist[which(dfmylist$refvec == "g"),"AA_category"]<-"Non-drastic, non-CpG" 
  dfmylist[which(dfmylist$refvec == "h"),"AA_category"]<-"Non-drastic, CpG" 
  dfmylist[which(dfmylist$refvec == "i"),"AA_category"]<-"Drastic, non-CpG" 
  dfmylist[which(dfmylist$refvec == "j"),"AA_category"]<-"Drastic, CpG" 
  
  dfmylist[which(dfmylist$refvec == "k"),"AA_category"]<-"Non-drastic, non-CpG" 
  dfmylist[which(dfmylist$refvec == "l"),"AA_category"]<-"Non-drastic, CpG" 
  dfmylist[which(dfmylist$refvec == "m"),"AA_category"]<-"Drastic, non-CpG" 
  dfmylist[which(dfmylist$refvec == "n"),"AA_category"]<-"Drastic, CpG"
  
  
  dfmylist[which(dfmylist$refvec == "o"),"AA_category"]<-"Non-drastic, non-CpG" 
  dfmylist[which(dfmylist$refvec == "p"),"AA_category"]<-"Drastic, non-CpG" 
  
  
  dfmylist[which(dfmylist$refvec == "q"),"AA_category"]<-"Non-drastic, non-CpG" 
  dfmylist[which(dfmylist$refvec == "r"),"AA_category"]<-"Drastic, non-CpG" 
  
  
  dfmylist[which(dfmylist$refvec == "a"),"color"]<-"green" 
  dfmylist[which(dfmylist$refvec == "b"),"color"]<-"blue" 
  dfmylist[which(dfmylist$refvec == "a"),"color"]<-"green" 
  dfmylist[which(dfmylist$refvec == "b"),"color"]<-"blue" 
  dfmylist[which(dfmylist$refvec == "c"),"color"]<-"green" 
  dfmylist[which(dfmylist$refvec == "d"),"color"]<-"blue" 
  dfmylist[which(dfmylist$refvec == "e"),"color"]<-"green" 
  dfmylist[which(dfmylist$refvec == "f"),"color"]<-"green" 
  
  dfmylist[which(dfmylist$refvec == "g"),"color"]<-"green" 
  dfmylist[which(dfmylist$refvec == "h"),"color"]<-"blue" 
  dfmylist[which(dfmylist$refvec == "i"),"color"]<-"orange" 
  dfmylist[which(dfmylist$refvec == "j"),"color"]<-"red" 
  
  dfmylist[which(dfmylist$refvec == "k"),"color"]<-"green" 
  dfmylist[which(dfmylist$refvec == "l"),"color"]<-"blue" 
  dfmylist[which(dfmylist$refvec == "m"),"color"]<-"orange" 
  dfmylist[which(dfmylist$refvec == "n"),"color"]<-"red"
  
  
  dfmylist[which(dfmylist$refvec == "o"),"color"]<-"green" 
  dfmylist[which(dfmylist$refvec == "p"),"color"]<-"orange" 
  
  
  dfmylist[which(dfmylist$refvec == "q"),"color"]<-"green" 
  dfmylist[which(dfmylist$refvec == "r"),"color"]<-"orange" 
  
  #View(dfmylist)
  #the plot will have a log y scale so all zero y values must go away.
  
  dfmylist <-dfmylist[which((dfmylist$values != 0)|(is.na(dfmylist$values))),]
  
  
  
  #split the data as there will be two plots one for "Synonomous" & one for "Non-synonomous" sites
  vec1<-c("a","b","c","d","e","f")
  dfmylist1<-dfmylist[which((dfmylist$refvec=="a")|(dfmylist$refvec=="b")|(dfmylist$refvec=="c")
                            |(dfmylist$refvec=="d")|(dfmylist$refvec=="e")|(dfmylist$refvec=="f")),]
  
  dfmylist2<-dfmylist[-which((dfmylist$refvec=="a")|(dfmylist$refvec=="b")|(dfmylist$refvec=="c")
                             |(dfmylist$refvec=="d")|(dfmylist$refvec=="e")|(dfmylist$refvec=="f")),]
  
  
  
  
  #prepare mean & error values for synonomous plot
  
  sem<-function(x){
    return(sd(x,na.rm = FALSE)/sqrt(length(x)))
  }
  
  dfmylist1$mean_vals = 0.0001
  dfmylist1$mean_vals[(which(dfmylist1$refvec == "a"))]<-mean(dfmylist1$values[(which(dfmylist1$refvec == "a"))], na.rm = FALSE)
  dfmylist1$mean_vals[(which(dfmylist1$refvec == "b"))]<-mean(dfmylist1$values[(which(dfmylist1$refvec == "b"))], na.rm = FALSE)
  dfmylist1$mean_vals[(which(dfmylist1$refvec == "c"))]<-mean(dfmylist1$values[(which(dfmylist1$refvec == "c"))], na.rm = FALSE)
  dfmylist1$mean_vals[(which(dfmylist1$refvec == "d"))]<-mean(dfmylist1$values[(which(dfmylist1$refvec == "d"))], na.rm = FALSE)
  dfmylist1$mean_vals[(which(dfmylist1$refvec == "e"))]<-mean(dfmylist1$values[(which(dfmylist1$refvec == "e"))], na.rm = FALSE)
  dfmylist1$mean_vals[(which(dfmylist1$refvec == "f"))]<-mean(dfmylist1$values[(which(dfmylist1$refvec == "f"))], na.rm = FALSE)
  
  dfmylist1$sem_vals = 0.0001
  dfmylist1$sem_vals[(which(dfmylist1$refvec == "a"))]<-sem(dfmylist1$values[(which(dfmylist1$refvec == "a"))])
  dfmylist1$sem_vals[(which(dfmylist1$refvec == "b"))]<-sem(dfmylist1$values[(which(dfmylist1$refvec == "b"))])
  dfmylist1$sem_vals[(which(dfmylist1$refvec == "c"))]<-sem(dfmylist1$values[(which(dfmylist1$refvec == "c"))])
  dfmylist1$sem_vals[(which(dfmylist1$refvec == "d"))]<-sem(dfmylist1$values[(which(dfmylist1$refvec == "d"))])
  dfmylist1$sem_vals[(which(dfmylist1$refvec == "e"))]<-sem(dfmylist1$values[(which(dfmylist1$refvec == "e"))])
  dfmylist1$sem_vals[(which(dfmylist1$refvec == "f"))]<-sem(dfmylist1$values[(which(dfmylist1$refvec == "f"))])
  
  sem(dfmylist1$values[(which(dfmylist1$refvec == "a"))])
  dfmylist1$UCLS = 0.0001
  dfmylist1$LCLS = 0.0001
  
  
  dfmylist1$LCLS = dfmylist1$mean_vals - dfmylist1$sem_vals
  dfmylist1$UCLS = dfmylist1$mean_vals + dfmylist1$sem_vals
  
  
  
  
  #create vector with named elements for ggplot to use to generate colors and legend
  col <- as.character(dfmylist1$color)
  col
  names(col) <- as.character(dfmylist1$AA_category)
  
  
  #plot synonomous points
  Synplot<-ggplot() +
    
    geom_point( data = dfmylist1, mapping = aes(x=refvec, y= values, colour = AA_category), size = 10,show.legend = FALSE) +
    xlab("Mutation Type") + ylab("Samples' Average Frequencies of Mutations") +
    geom_point(data = dfmylist1, mapping = aes (x = refvec, y = mean_vals, colour = AA_category),size = 5.0) +
    geom_errorbar(data = dfmylist1,aes(x = refvec, ymin= LCLS, ymax= UCLS, color = AA_category),width=.5) +
    scale_color_manual(values=col) +
    theme(legend.text=element_text(size=7),legend.title = element_text(size=10)) +
    
    
    
    
    annotate("text", x  = 1.5, y = 0.00001, label = "A -> G") +
    
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(2.5)) +
    
    annotate("text", x  = 3.5, y = 0.00001, label = "T -> C") +
    
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(4.5)) +
    
    geom_point(data = dfmylist1, mapping = aes (x = "e1", y = 0.0), colour = "red",size = 0.0) +
    
    geom_point(data = dfmylist1, mapping = aes (x = "f1", y = 0.0), colour = "red",size = 0.0) +
    
    annotate("text", x  = 5.5, y = 0.00001, label = "C -> T") +
    
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(6.5)) +
    
    annotate("text", x  = 7.5, y = 0.00001, label = "G -> A") +
    
    scale_x_discrete(labels=c("a" = "", "b" = "", "c" = "", "d" = "","e" = "", "e1" = "","f" ="", "f1" = "")) +
    
    scale_y_log10(labels = comma) +
    theme(legend.position="none") +
    expand_limits(y = c(0.00001, 0.1)) +
    theme(plot.margin = unit(c(1,1,3.0,1), "cm")) +
    theme (panel.border = element_rect(colour = "black", fill=NA, size=3),plot.title = element_text(hjust = 0.5))
  
  
  
  
  
  
  #prepare mean & error values for nonsynonomous plot
  
  dfmylist2$mean_vals = 0.0001
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "g"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "g"))], na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "h"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "h"))],na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "i"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "i"))], na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "j"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "j"))], na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "k"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "k"))], na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "l"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "l"))], na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "m"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "m"))], na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "n"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "n"))], na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "o"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "o"))], na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "p"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "p"))], na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "q"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "q"))], na.rm = FALSE)
  dfmylist2$mean_vals[(which(dfmylist2$refvec == "r"))]<-mean(dfmylist2$values[(which(dfmylist2$refvec == "r"))], na.rm = FALSE)
  
  
  dfmylist2$sem_vals = 0.0001
  
  
  
  sem(dfmylist2$values) 
  sem(dfmylist2$values[(which(dfmylist2$refvec == "g"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "g"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "g"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "h"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "h"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "i"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "i"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "j"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "j"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "k"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "k"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "l"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "l"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "m"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "m"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "n"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "n"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "o"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "o"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "p"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "p"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "q"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "q"))])
  dfmylist2$sem_vals[(which(dfmylist2$refvec == "r"))]<-sem(dfmylist2$values[(which(dfmylist2$refvec == "r"))])
  
  dfmylist2$UCLS = 0.0001
  dfmylist2$LCLS = 0.0001
  
  
  dfmylist2$LCLS = dfmylist2$mean_vals -dfmylist2$sem_vals
  dfmylist2$UCLS = dfmylist2$mean_vals + dfmylist2$sem_vals
  
  
  
  
  
  #create vector with named elements for ggplot to use to generate colors and legend
  col <- as.character(dfmylist2$color)
  col
  names(col) <- as.character(dfmylist2$AA_category)
  
  #plot synonomous points
  NonSynplt<-ggplot() +
    
    #plot points
    geom_point( data = dfmylist2, mapping = aes(x=refvec, y= values, colour = AA_category), size = 0.1,show.legend = TRUE) +
    xlab("Mutation Type") + ylab("Samples' Average Frequencies of Mutations") +
    geom_point(data = dfmylist2, mapping = aes (x = refvec, y = mean_vals, colour = AA_category),size = 5.0) +
    geom_errorbar(data = dfmylist2,aes(x = refvec, ymin= LCLS, ymax= UCLS, color = AA_category),width=.5) +
    scale_color_manual(values=col) +
    theme(legend.text=element_text(size=7),legend.title = element_text(size=10)) +
    
    #set up graph background
    annotate("text", x  = 2.5, y = 0.00001, label = "A -> G") +
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(4.5)) +
    annotate("text", x  = 6.5, y = 0.00001, label = "T -> C") +
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(8.5)) +
    annotate("text", x  = 9.5, y = 0.00001, label = "C -> T") +
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(10.5)) +
    annotate("text", x  = 11.5, y = 0.00001, label = "G -> A") +
    scale_x_discrete(labels=c("g" = "", "h" = "", "i" = "", "j" = "","k" = "","l" ="","m" ="", "n"="","o" = "","p" = "","q"="","r"="")) +
    
    #scales, margins, border
    scale_y_log10(labels = comma) +
    expand_limits(y = c(0.00001, 0.1)) +
    theme(plot.margin = unit(c(1,1,3.0,1), "cm")) +
    theme (panel.border = element_rect(colour = "black", fill=NA, size=3),plot.title = element_text(hjust = 0.5))
  
  
  
 # prints both graphs seperately 
  
  
  print(NonSynplt)
  
 print(Synplot)


}


Fig3(dfx)

