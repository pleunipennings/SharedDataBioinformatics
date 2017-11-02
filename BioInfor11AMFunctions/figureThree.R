figureThree <-function(dfx){
  
  CSV<-dfx
  
  class(CSV)
  
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
  mylist <- list (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) 
  mylist
  
  # Finding maximum row length and filling empty space with "NA"
  llngths<-lapply(mylist, function(x) length(x))
  vlngths<-unlist(llngths)
  maxlngth <- max(vlngths)
  
  myfn<-function(x,maxlngth){
    length(x)<-maxlngth
    return (x)}
  mtxmylst <-sapply(mylist, myfn, maxlngth)
  mtxmylst
  
  rm (mxstats)
  fnclcls<-function(x){
    if(mean(x, na.rm = TRUE) - 0.3*sd(x, na.rm = TRUE)<0){
      return (mean(x, na.rm = TRUE) - 0.001)}
    return(mean(x, na.rm = TRUE) - 0.3*sd(x, na.rm = TRUE))}
  means<-sapply(mylist, function(x) mean(x,na.rm = TRUE))
  ucls<-sapply(mylist, function(x) (0.3*sd(x, na.rm = TRUE) + mean(x, na.rm = TRUE)))
  lcls <- sapply (mylist, fnclcls)
  mxstats<-rbind(means,ucls)
  mxstats<-rbind(mxstats,lcls)
  mxstats
  
  # Constructing the data frame
  dfval<-data.frame(mtxmylst)
  names<-c("means","ucls","lcls")
  dfstats<-data.frame(mxstats,row.names = names)
  dfstats
  
  # Labelling/Naming
  rename
  dfval<-rename(dfval,c("X1"="val_a","X2"="val_b","X3"="val_c","X4"="val_d","X5"="val_e","X6"="val_f","X7"="val_g","X8"="val_h","X9"="val_i","X10"="val_j","X11"="val_k","X12"="val_l","X13"="val_m","X14"="val_n","X15"="val_o","X16"="val_p","X17"="val_q", "X18"= "val_r"))
  dfstats<-rename(dfstats,c("X1"="mut_a","X2"="mut_b","X3"="mut_c","X4"="mut_d","X5"="mut_e","X6"="mut_f","X7"="mut_g","X8"="mut_h","X9"="mut_i","X10"="mut_j","X11"="mut_k","X12"="mut_l","X13"="mut_m","X14"="mut_n","X15"="mut_o","X16"="mut_p","X17"="mut_q", "X18"= "mut_r"))
  dfval[,"val_m"]
  mean(dfval[,"val_m"])
  dfstats
  sink("Routput(sap).txt")
  print(dfstats)
  sink()
  
  class(dfval[,2])
  class(dfstats[,2])
  
  # Creating Synonymous Plot
  Synplt <- ggplot() +
    
    geom_point( data = dfval, mapping = aes(x="a", y= val_a), colour = "green", size = 0.1) +xlab(" Basepair") + ylab(" Frequencies of Mutations") +
    geom_point(data = dfstats, mapping = aes (x = "a", y = dfstats["means","mut_a"]), colour = "green",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "a", ymin= dfstats["lcls","mut_a"], ymax= dfstats["ucls","mut_a"]), color = "green",width=.5) +
    
    annotate("text", x  = 1.5, y = 0.00001, label = "A -> G") +
    
    geom_point( data = dfval, mapping = aes(x="b", y= val_b), colour = "blue", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "b", y = dfstats["means","mut_b"]), colour = "blue",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "b", ymin= dfstats["lcls","mut_b"], ymax= dfstats["ucls","mut_b"]), color = "blue",width=.5) +
    
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(2.5)) +
    
    geom_point( data = dfval, mapping = aes(x="c", y= val_c), colour = "green", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "c", y = dfstats["means","mut_c"]), colour = "green",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "c", ymin= dfstats["lcls","mut_c"], ymax= dfstats["ucls","mut_c"]), color = "green",width=.5) +
    
    annotate("text", x  = 3.5, y = 0.00001, label = "T -> C") +
    
    geom_point( data = dfval, mapping = aes(x="d", y= val_d), colour = "blue", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "d", y = dfstats["means","mut_d"]), colour = "blue",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "d", ymin= dfstats["lcls","mut_d"], ymax= dfstats["ucls","mut_d"]), color = "blue",width=.5) +
    
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(4.5)) +
    
    geom_point( data = dfval, mapping = aes(x="e", y= val_e), colour = "green", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "e", y = dfstats["means","mut_e"]), colour = "green",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "e", ymin= dfstats["lcls","mut_e"], ymax= dfstats["ucls","mut_e"]), color = "green",width=.5) +
    
    geom_point(data = dfstats, mapping = aes (x = "e1", y = 0.0), colour = "red",size = 0.0) +
    
    annotate("text", x  = 5.5, y = 0.00001, label = "C -> T") +
    
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(6.5)) +
    
    geom_point( data = dfval, mapping = aes(x="f", y= val_f), colour = "green", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "f", y = dfstats["means","mut_f"]), colour = "green",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "f", ymin= dfstats["lcls","mut_f"], ymax= dfstats["ucls","mut_f"]), color = "green",width=.5) +
    
    geom_point(data = dfstats, mapping = aes (x = "f1", y = 0.0), colour = "red",size = 0.0) +
    
    annotate("text", x  = 7.5, y = 0.00001, label = "G -> A") +
    
    scale_x_discrete(labels=c("a" = "", "b" = "", "c" = "", "d" = "","e" = "", "e1" = "","f" ="", "f1" = "")) +
    
    #ggtitle("HIV Genetic Mutations Study - Frequency vs Type") +
    scale_y_log10(labels = comma) +
    theme(legend.position="none") +
    expand_limits(y = c(0.00001, 0.1)) +
    theme(plot.margin = unit(c(1,1,3.0,1), "cm")) +
    theme (panel.border = element_rect(colour = "black", fill=NA, size=3),plot.title = element_text(hjust = 0.5))
  
  # Creating nonsynonymous Plot
  NonSynplt <- ggplot() +
    
    geom_point( data = dfval, mapping = aes(x="g", y= val_g), colour = "green", size = 0.1) +xlab(" Basepair") + ylab(" Frequencies of Mutations") +
    geom_point(data = dfstats, mapping = aes (x = "g", y = dfstats["means","mut_g"]), colour = "green",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "g", ymin= dfstats["lcls","mut_g"], ymax= dfstats["ucls","mut_g"]), color = "green",width=.5) +
    
    annotate("text", x  = 2.5, y = 0.00001, label = "A -> G") +
    
    geom_point( data = dfval, mapping = aes(x="h", y= val_h), colour = "blue", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "h", y = dfstats["means","mut_h"]), colour = "blue",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "h", ymin= dfstats["lcls","mut_h"], ymax= dfstats["ucls","mut_h"]), color = "blue",width=.5) +
    
    
    geom_point( data = dfval, mapping = aes(x="i", y= val_i), colour = "orange", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "i", y = dfstats["means","mut_i"]), colour = "orange",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "i", ymin= dfstats["lcls","mut_i"], ymax= dfstats["ucls","mut_i"]), color = "orange",width=.5) +
    
    geom_point( data = dfval, mapping = aes(x="j", y= val_j), colour = "red", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "j", y = dfstats["means","mut_j"]), colour = "red",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "j", ymin= dfstats["lcls","mut_j"], ymax= dfstats["ucls","mut_j"]), color = "red",width=.5) +
    
    
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(4.5)) +
    
    
    annotate("text", x  = 6.5, y = 0.00001, label = "T -> C") +
    
    geom_point( data = dfval, mapping = aes(x="k", y= val_k), colour = "green", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "k", y = dfstats["means","mut_k"]), colour = "green",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "k", ymin= dfstats["lcls","mut_k"], ymax= dfstats["ucls","mut_k"]), color = "green",width=.5) +
    
    
    geom_point( data = dfval, mapping = aes(x="l", y= val_l), colour = "blue", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "l", y = dfstats["means","mut_l"]), colour = "blue",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "l", ymin= dfstats["lcls","mut_l"], ymax= dfstats["ucls","mut_l"]), color = "blue",width=.5) +
    
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(8.5)) +
    
    geom_point( data = dfval, mapping = aes(x="m", y= val_m), colour = "orange", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "m", y = dfstats["means","mut_m"]), colour = "orange",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "m", ymin= dfstats["lcls","mut_m"], ymax= dfstats["ucls","mut_m"]), color = "orange",width=.5) +
    
    annotate("text", x  = 9.5, y = 0.00001, label = "C -> T") +
    
    geom_point( data = dfval, mapping = aes(x="n", y= val_n), colour = "red", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "n", y = dfstats["means","mut_n"]), colour = "red",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "n", ymin= dfstats["lcls","mut_n"], ymax= dfstats["ucls","mut_n"]), color = "red",width=.5) +
    
    geom_vline(aes(linetype=1, colour="black"),xintercept =c(10.5)) +
    
    geom_point( data = dfval, mapping = aes(x="o", y= val_o), colour = "green", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "o", y = dfstats["means","mut_o"]), colour = "green",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "o", ymin= dfstats["lcls","mut_o"], ymax= dfstats["ucls","mut_o"]), color = "green",width=.5) +
    
    annotate("text", x  = 11.5, y = 0.00001, label = "G -> A") +
    
    
    geom_point( data = dfval, mapping = aes(x="p", y= val_p), colour = "orange", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "p", y = dfstats["means","mut_p"]), colour = "orange",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "p", ymin= dfstats["lcls","mut_p"], ymax= dfstats["ucls","mut_p"]), color = "orange",width=.5) +
    
    
    geom_point( data = dfval, mapping = aes(x="q", y= val_q), colour = "green", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "q", y = dfstats["means","mut_q"]), colour = "green",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "q", ymin= dfstats["lcls","mut_q"], ymax= dfstats["ucls","mut_q"]), color = "green",width=.5) +
    
    
    geom_point( data = dfval, mapping = aes(x="r", y= val_r), colour = "orange", size = 0.1) +
    geom_point(data = dfstats, mapping = aes (x = "r", y = dfstats["means","mut_r"]), colour = "orange",size = 5.0) +
    geom_errorbar(data = dfstats,aes(x = "r", ymin= dfstats["lcls","mut_r"], ymax= dfstats["ucls","mut_r"]), color = "orange",width=.5) +
    
    scale_x_discrete(labels=c("g" = "", "h" = "", "i" = "", "j" = "","k" = "","l" ="","m" ="", "n"="","o" = "","p" = "","q"="","r"="")) +
    
    #ggtitle("HIV Genetic Mutations Study - Frequency vs Type") +
    scale_y_log10(labels = comma) +
    theme(legend.position="none") +
    expand_limits(y = c(0.00001, 0.1)) +
    theme(plot.margin = unit(c(1,1,3.0,1), "cm")) +
    theme (panel.border = element_rect(colour = "black", fill=NA, size=3),plot.title = element_text(hjust = 0.5))
  
  #, vp=viewport(width=1.0, height=0.97)
  require(grid)
  require(gridExtra)
  title1=textGrob("
                  Fig. 3: type of site vs frequency", gp=gpar(fontface="bold", fontsize = 16, cex = 1))
  grid.arrange( top =title1,Synplt + ggtitle('Synonymous Plots'), NonSynplt + ggtitle('nonsynonymous Plots'),  nrow=1)
  
  # Creating the plot keys
  grid.text("KEY:
            Color Green = Non Drastic AA change, non-CpG forming
            Color Blue = Non Drastic AA change, CpG forming
            Color Orange = Drastic AA change, non-CpG-forming
            Color Red = Drastic AA change, CpG forming ", 
            x = unit(2, "cm"), y = unit(0.25,"cm"), just = "left", vjust = unit(0.0,"cm"))
}
figureThree(yourDF)
