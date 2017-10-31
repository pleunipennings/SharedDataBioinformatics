#Nathan
#### in class coursework stuff! ####
h1 <- read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05",format="fasta")
h1
#?read.alignment
h2<- dist.alignment(h1)
h2
nj(h1)
?nj
h3<-nj(h2
)
plot(h3)


#nate stuff
G <- read.alignment("class25Influ.txt", format = "fasta")
G

G1 <- dist.alignment(G)
G1

G2 <- nj(G1)
G2

G3 <- plot(G2)
G3

?`ape-package`
library(help = "ape")


ahub <- AnnotationHub()
ahub["AH5086"]

#### group notes ####


seq_GC <- read.fasta("class25Influ.txt", format = "fasta")
dist.alignment(seq_GC)
x= split(seq, (0:nrow(seq) %/% 500))

#?gregexpr
seqT = data.frame (x = "tttttattttccaactactggcattggcatcggcatccataattaccttt
                   gtgaaacttgcacaatccttgttctacctttgtttttttttctttttaca
                   aaatgtacttagtattagaagactatgactatctgatctaaaccaaactt
                   catagacatgtgtggtttattagcttgaatgctaacgtgcagatgaatga
                   agtggtagaggaaaaatgacagaattgctttcctccttttttatcataca
                   tctgtggaaaaatcctgaagtgttacacagctctatcatgtcatcaaata
                   gtttatgtataactctactattgtctgagggtctgtgtgcattgatgata
                   tgcccatgaaacaggcacagtttctttttctttttttctgtttttttttt
                   tttttttttttttttgagatggagtcttgctctgttgcccaggctggagt
                   gcagtggagtggcgcaatatgggctcactgcaacctccacctcccaggtt
                   caagtgattctcctgcgtcagcctcccaggtagctgggattacaggcatg
                   ccccaccatgcccggctaatttttgttttttgtttgtttgtttgtttgtt
                   ttgtttttttttgagacagagtctcgctctgtcaccctggctggagtgca
                   gtggcgcgttctaggctcactgcagactccacctcccgggttcacgccat
                   tctcctggctcagcctccggagtagctgggactacaggcacccaccacca
                   cgcctggctaattttgtttttgtaattttagtagaggcggggtttcaccg
                   tgttagccaggatggtctcgatctcctgaccttgtgatccgccctcctcg
                   gcctcccaaagtgctgggattacaggcatgagccaccacgcccggcctaa
                   tttttgtattttttaatagagacggggtttcaccatgttggtcaggctgg
                   tctcgaactcctgacctcaggtgattcacccgcctcagcctcccaaagtg
                   ctgggattacaggcgtgagccaatgcgcccggcccaacaggcacagtttc
                   tgctctgatggaacttacagtctaacaagggataaccaggatgttcatga
                   aatttgtagacataagcaaatatatgtataatgtattaagaatgtctaca
                   gttggttaaaaaaataaagtgatcttgtctctatttccaagaccatttca
                   atgtcatatatttagttgtgcaacagtgggacaatggggtactatgtgaa
                   cagctaactgtgtggaatgagattcaaagaatgtttgtttgaagaggtga
                   tgtctaagctaggtaactctagaagagtaagcagctctctctctaggaag
                   agcgttctagtcagaaagaagagcatgtgctgaagggcccagagaaagca
                   tgacacattcaggaactgaaaatatgtcagtatggctaaagcacagaatt
                   ccaggggaaagatctattgtgaaataataggctggagagagaagcaggaa
                   gcagattgtgcaaggtttttaaaacccaagttaaggaatttagaggtcag
                   gggatttggcctgggggtgacataatcagattgtcgcttagaaagatgat
                   tcttgctgtggaacagagtatggattggtgcgatgtgggagtgacaggag
                   aaaaacctgttgcagtaatctaaaccagcaatgtggtggcttgaattagc
                   agagtagtaagggagataagtggagagaatttagaaatatttaggaagtg
                   gatttgacaggcctcagtgattcatttgatgtggtgcagggatgtccaat
                   cttttggcttctccgggccacagtggaagaagaattgtcttggaccccac
                   ataaaatacactaacgctaacggtagctgatgagctaaaaaaaaaaaaaa
                   aaaaaaaaaaatcgcaaaaaaaaaaaaaaatgttttaggaaagtttacga
                   attcatttgggcctcattcaaagccatcctgggccacaggttggacaagc
                   ttgatgtggtgggtgaaagggagtcaaggaagacacttaggtgcttggct
                   tggaatatcaggtgaatgtgggacttgtctggtaaggttagacacggagc
                   tgtggtgctcaggagaaatatacaaacaacgtttagggatcatgagcctg
                   taggtgctgcttgtaacgggggagaggattggatagctccaggaaagtgt")

a=lapply(x, function(vec){
    x <- gregexpr("gc", vec, perl = TRUE)
    res <- sum(attr(x[[1]], "match.length"))
    res
})

b=lapply(x, function(vec){
    x <- gregexpr("g", vec, perl = TRUE)
    res <- regmatches(vec, x)
    res
}) 


c=lapply(x, function(vec){
    x <- gregexpr("c", vec, perl = TRUE)
    res <- regmatches(vec, x)
    res
}) 

slide_function <- function(data, window, step){
    total <- length(data)
    spots <- seq(from = 1, to = (total-window), length.out = step)
    result <- vector(length = length(spots))
    for(i in 1:length(spots)){result[i]}
    return(result)
}
slide_function(seqT,2,500)
?seq.default
#### from divergence class data ####
as.table(HPIV1a, stringsAsFactors=FALSE)
M <- as.data.frame(read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta"))
# Here's a function that we've coded that will calculate divergence for basepair windows 
#   across the chromosome. We create it using a function called "function", and assign it to 
#   it's own function named "div.by.window".
#   The output is a dataframe with 5 columns, 
#     1) the number of the window (in sequential order)
#     2) the starting chromosomal position
#     3) the ending chromosomal position
#     4) the number of bps in the window
#     5) the proportion of diverged sites in the window
#   Run this code in order to use the function.
div.by.window = function(x, window_size){
    # window_size = number of bp per window
    # x = dataframe, with a column, "div" with a TRUE value for a diverged site and FALSE for non-diverged
    window_start = seq(from=1, to=nrow(x), by=window_size)
    diverged.window = data.frame(window=1:length(window_start), start.pos=NA, end.pos=NA, length=NA, p.diverged=NA)
    for (i in 1:length(window_start) ){
        focal_start = window_start[i]
        if (i == length(window_start) ){
            focal_end = nrow(x)
        } else {
            focal_end = focal_start + window_size  
        }
        focal = x[focal_start:focal_end, ]
        diverged.window[i,2] = focal[1,1]
        diverged.window[i,3] = focal[nrow(focal),1]
        diverged.window[i,4] = diverged.window[i,3]-diverged.window[i,2]
        if ( nrow(focal) <= 0 ) {next}
        diverged.window[i,5] = mean(focal$div)  
    }
    diverged.window
}

# Use this new function "div.by.window". It takes two arguments, "x" (your humanchimp data object)
#   and "window_size" (the number of bp that you want per window). For now, use 5000 bp, 
#   which is the approximate size of a small gene. 
#   Assign the output to the object "diverged.window".
diverged.window <- div.by.window(x = humanchimp,
                                 window_size = 5000
)
mean(diverged.window$p.diverged)
#### python not R code does not work for us ####
# below is a non R code for a CPG Island Finder with Sliding Window Algorithm: 
    #String Index Out of Bound Exception intermittently

public static List<Integer> finalCPGIslands(List<Integer> iList,
                                            String iSeq, int width) {
    # Declare output list that contains final list of start and end
    # intervals
    List<Integer> oList = new ArrayList<Integer>();
    # Add the first two elements anyways
    oList.add(iList.get(0));
    oList.add(iList.get(1));
    if (iList.size() > 2) {
        for (int i = 2; i < iList.size(); i += 2) {
            # The below IF is attempted to ensure that substring is always
            # valid
            if (iSeq.length() > iList.get(i + 1)) {
                # While creating the substring in next line, I get String
                # index out of range: -9
                String testSeq = iSeq.substring(iList.get(i),
                                                iList.get(i + 1) + 1);
                boolean check = cpgCriteriaCheck(testSeq);
                if (check) {
                    # If condition is met, add the indexes to the final
                    #list
                    oList.add(iList.get(i));
                    oList.add(iList.get(i + 1));
                }
                # If condition is not met, start removing one character at
                # a time until condition is met
                else {
                    
                    int counter = 0;
                    int currentSequenceLength = testSeq.length();
                    String newTestSeq = null;
                    while (counter <= currentSequenceLength) {
                        counter++;
                        if (testSeq.length() > 2) {
                            newTestSeq = testSeq.substring(1,
                                                           testSeq.length() - 1);
                            testSeq = newTestSeq;
                            if (newTestSeq.length() < width) {
                                counter = currentSequenceLength + 1;
                            } else {boolean checkAgain = cpgCriteriaCheck(newTestSeq);
                            # If condition met, add the item to list
                            # and exit
                            if (checkAgain) {
                                oList.add(iList.get(i) + counter);
                                oList.add(iList.get(i + 1) - counter);
                                counter = currentSequenceLength + 1;
                            }
                            
                            } # End of Else
                        } # End of IF
                        
                    } # End of While
                } # End of Else
            }
            
        } # End of For
    } # End of Else
    return oList;
    #example("cpg.assoc")
    #biostrings
    
    library("Biostrings")
    
    s = readDNAStringSet("nm.fasta")
    subseq(s, start=c(1, 2, 3), end=c(3, 6, 5))
    
    ## try http:// if https:// URLs are not supported
    # source("https://bioconductor.org/biocLite.R")
    #biocLite("AnnotationHub")
    
    

################## freq finder? maybe? ####################
    
    import sys,re,fileinput
    
    Argument = []
    Argument = sys.argv[1:] 
    
    if (len(Argument)) < 3:
      print "Usage: Input_pileup_file Column_with_information_about_read_bases(usually_column_9_in_pileup) Output_file" 
    sys.exit()
    
    File_Pileup = Argument[0]
    index = int(Argument[1])-1
    output = open(str(Argument[2]),"w")
    
    
    nucleotides = ["A","T","C","G","a","t","c","g"]
    complement = {'a':'T','g':'C','t':'A','c':'G'}
    
    output.write("Chromosome\tCoordinate\tReference_base")
    
    Frequency_bases = {"A":0,"T":0,"C":0,"G":0}
    
    for base in sorted(Frequency_bases.keys()):
      output.write("\t"+str(base))
    
    output.write("\n")
    
    for line in fileinput.input([File_Pileup]):
      rowlist = []
    rowlist = (line.rstrip("\n")).split('\t')
    rowlist[2] = rowlist[2].upper()
    
    Frequency = {"A":0,"T":0,"C":0,"G":0}
    
    if line.startswith("#"):
      continue
    else:
      output.write(str(rowlist[0])+"\t"+str(rowlist[1])+"\t"+str(rowlist[2]))
    for i in rowlist[index]:
      if i == ".":
      Frequency[rowlist[2]] = Frequency[rowlist[2]] + 1
    continue
    if i == ",":
      Frequency[rowlist[2]] = Frequency[rowlist[2]] + 1
    continue
    
    if i in nucleotides:
      if i.isupper():
      Frequency[i] = Frequency[i]+1
    
    else:
      Frequency[complement[i]] = Frequency[complement[i]] + 1
    
    for base in sorted(Frequency.keys()):
      output.write("\t"+str(Frequency[base]))
    
    output.write("\n")
    
    output.close()
#### Mergerger ####
    
    HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
    consensus(HPIV1a)
    Error in obj[[1]]$tip.label : $ operator is invalid for atomic vectors
    ??read.alignment
    ??consensus
    example(ape::consensus)
    seqinr::consensus(HPIV1a)
    
length(seqinr::consensus(HPIV1a))
    [1] 15473
    

?plot.default


plot.default(x = c(V1, V2), #plot it!
             xlab = "Location", ylab = "Frequency", main = "Location vs frequency CpG non-CpG Graph",
             col =  c("blue","red")
)
read.csv("OverviewSelCoeff_BachelerFilter.csv")->n
head(n)

#follows is reinterpretation of Pleuni's fig2 code
#im trying to figure out how the plot is formatted....
LvsF_CpG_Printer <- function(data_frame){
  png("Location_vs_frequrency_CpG_non-CpG_Graph.png",width=12,height=7.5,units="in",res=100)
  par(mfrow=c(1,1))
  maxnuc=1000
  par(mar = c(3,5,1,2))
  plot(n$num[40:maxnuc],n$EstSelCoeff[40:maxnuc],
       log="y", ylab="Frequrency",cex.lab=1.3,
       xaxt="n",yaxt="n", xlab="",
       col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-4,1),xlim=c(40,979))
        axis(1,at=c(3*seq(15,95,by=20)-1,296+30),labels=c(seq(15,95,by=20),""))
        axis(1,at=3*seq(109,349,by=20)-1,labels=seq(109-99,349-99,by=20))
        axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),las=1,line=0,tick=FALSE)
        eaxis(side = 2, at = 10^((-0):(-(5))),label=rep("",6))
         rect(0, 0.00001, 297.5, 2, density = NULL, angle = 45,col="grey70",border = NA)
          for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray41")}
         cols <- brewer.pal(6, "Set2")[c(1, 2, 3, 6)]
          for (i in 40:maxnuc){
            c=0; co = 1
            if (n$TypeOfSite[i]=="stop"&n$WTnt[i]%in%c("g","c")) {c=1;p=21}
            if (n$TypeOfSite[i]=="syn"&n$WTnt[i]%in%c("g","c")) {c=cols[1];p=21}
            if (n$TypeOfSite[i]=="syn"&n$WTnt[i]%in%c("a","t")) {c=cols[1];p=21}
            if (n$TypeOfSite[i]=="nonsyn"&n$WTnt[i]%in%c("c","g")) {c=cols[2];p=21}
            if (n$TypeOfSite[i]=="nonsyn"&n$WTnt[i]%in%c("a","t")) {c=cols[4];p=21}
            if (i %in% 73:81) {p = 22; co = 2} 
            if (c!=0) points(n$num[i],n$EstSelCoeff[i],pch=p,col=co,
                     bg=rgb(red=col2rgb(c)[1]/255,
                            green=col2rgb(c)[2]/255,
                            blue=col2rgb(c)[3]/255,
                            maxColorValue = 1,alpha=0.8),
                     cex=2)
            }
          rect(0, 0.000001, 1200, 3.5*10^-4, density = NULL, angle = 45,col=1,border = NA)
          text(55*3,2.9*10^-4,"PROTEASE",col="white")
          rect(297.5, 0.000001, 1200, 3.5*10^-4, density = NULL, angle = 45,col="grey40",border = NA)
          text(220*3,2.9*10^-4,"REVERSE TRANSCRIPTASE",col="white")
            legpos=296;legposV=0.4
            rect(legpos*3, 0.4*legposV, (legpos+42.5)*3, 1.7*legposV, density = NULL, angle = 45,col=alpha("white",1))
            points((legpos+5)*3,legposV/0.7,pch=21,bg=1,col=1,cex=2)
            text((legpos+9)*3,legposV/0.7,"Nonsense",adj=0)
            points((legpos+5)*3,legposV,pch=21,bg=cols[2],col=1,cex=2)
            text((legpos+9)*3,legposV,"Non-syn, C/G",adj=0)
            points((legpos+5)*3,legposV*0.7,pch=21,bg=cols[4],col=1,cex=2)
            text((legpos+9)*3,legposV*0.7,"Non-syn, A/T",adj=0)
            points((legpos+5)*3,legposV*0.49,pch=21,bg=cols[1],col=1,cex=2)
            text((legpos+9)*3,legposV*0.49,"Synonymous",adj=0)
dev.off()
return()        
}

LvsF_CpG_Printer(n)




FreqsSyn<-n$MeanFreq[n$TypeOfSite=="syn"]
FreqsNonSyn<-n$MeanFreq[n$TypeOfSite=="nonsyn"]
FreqsStop<-n$MeanFreq[n$TypeOfSite=="stop"]

wilcox.test(FreqsSyn, FreqsNonSyn,alternative = "greater", paired = FALSE)
wilcox.test(FreqsNonSyn,FreqsStop,alternative = "greater", paired = FALSE)

if (TRUE){
  #pdf("../Output/EstSelCoeffPRO_aug2017.pdf",width=12,height=8)
  png("Location_vs_frequrency_CpG_non-CpG_Graph.png",width=12,height=7.5,units="in",res=100)
  #?png #prints out the below plot as a png!
  par(mfrow=c(1,1))
  maxnuc=1000
  #??maxnuc #maxnuc? maxnuc sets the maximum number of variables to be used in 
  # the plot
  par(mar = c(3,5,1,2))
  #?par #what does par do?
  plot(n$num[40:maxnuc],n$EstSelCoeff[40:maxnuc],
       log="y", ylab="Frequrency",cex.lab=1.3,
       xaxt="n",yaxt="n", xlab="",
       col="darkgrey",t="n",pch=".", ylim=c(3.2*10^-4,1),xlim=c(40,979))
  
  axis(1,at=c(3*seq(15,95,by=20)-1,296+30),labels=c(seq(15,95,by=20),""))
  axis(1,at=3*seq(109,349,by=20)-1,labels=seq(109-99,349-99,by=20))
  axis(2,at=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),labels=c(10^-5,10^-4,10^-3,10^-2,10^-1,10^-0),las=1,line=0,tick=FALSE)
  eaxis(side = 2, at = 10^((-0):(-(5))),label=rep("",6))
  
  #color Protease region grey
  rect(0, 0.00001, 297.5, 2, density = NULL, angle = 45,col="grey70",border = NA)
  for(i in 1:5){abline(h = 1:10 * 10^(-i), col = "gray41")}
  
  cols <- brewer.pal(6, "Set2")[c(1, 2, 3, 6)]
  for (i in 40:maxnuc){
    c=0; co = 1
    if (n$TypeOfSite[i]=="stop"&n$WTnt[i]%in%c("g","c")) {c=1;p=21}
    if (n$TypeOfSite[i]=="syn"&n$WTnt[i]%in%c("g","c")) {c=cols[1];p=21}
    if (n$TypeOfSite[i]=="syn"&n$WTnt[i]%in%c("a","t")) {c=cols[1];p=21}
    if (n$TypeOfSite[i]=="nonsyn"&n$WTnt[i]%in%c("c","g")) {c=cols[2];p=21}
    if (n$TypeOfSite[i]=="nonsyn"&n$WTnt[i]%in%c("a","t")) {c=cols[4];p=21}
    if (i %in% 73:81) {p = 22; co = 2} #for Active site Protease change pch
    if (c!=0) points(n$num[i],n$EstSelCoeff[i],pch=p,col=co,
                     bg=rgb(red=col2rgb(c)[1]/255,
                            green=col2rgb(c)[2]/255,
                            blue=col2rgb(c)[3]/255,
                            maxColorValue = 1,alpha=0.8),
                     cex=2)
  }
  
  #Add "Protease" and "RT" words
  rect(0, 0.000001, 1200, 3.5*10^-4, density = NULL, angle = 45,col=1,border = NA)
  text(55*3,2.9*10^-4,"PROTEASE",col="white")
  rect(297.5, 0.000001, 1200, 3.5*10^-4, density = NULL, angle = 45,col="grey40",border = NA)
  text(220*3,2.9*10^-4,"REVERSE TRANSCRIPTASE",col="white")
  
  
  
  #Add legend
  legpos=296;legposV=0.4
  rect(legpos*3, 0.4*legposV, (legpos+42.5)*3, 1.7*legposV, density = NULL, angle = 45,col=alpha("white",1))
  points((legpos+5)*3,legposV/0.7,pch=21,bg=1,col=1,cex=2)
  text((legpos+9)*3,legposV/0.7,"Nonsense",adj=0)
  points((legpos+5)*3,legposV,pch=21,bg=cols[2],col=1,cex=2)
  text((legpos+9)*3,legposV,"Non-syn, C/G",adj=0)
  points((legpos+5)*3,legposV*0.7,pch=21,bg=cols[4],col=1,cex=2)
  text((legpos+9)*3,legposV*0.7,"Non-syn, A/T",adj=0)
  points((legpos+5)*3,legposV*0.49,pch=21,bg=cols[1],col=1,cex=2)
  text((legpos+9)*3,legposV*0.49,"Synonymous",adj=0)
  
  ?dev.off()
}

  #  1/ make consensus data of our sample DNA (for one file only!)
n <- data.frame(Pos = c(1:length(seqinr::consensus(HPIV1a))),
                WTnt = (seqinr::consensus(HPIV1a)),
                Trans = c(1)
                )
rho((seqinr::consensus(HPIV1a)), wordsize = 2, alphabet = s2c("cg"))

head(n)
tail(n)
n$Trans
?seqinr::consensus
alphabetFrequency(n$WTnt, c("c","g","a","t"), as.prob = T, baseOnly=TRUE
)
consensusString(HPIV1a, baseOnly=TRUE, ambiguityMap=IUPAC_CODE_MAP)

  #  2/ determine the frequency of mutation at each base
  #  3/ create a function that will plot the cpg /posible cpg

zscore((seqinr::consensus(HPIV1a)), modele = "codon")




################### pruned section from master script ###################
#setwd(SharedDataBioinformatics/11AM_Influ_team)    
seq_CG <- read.alignment("class25Influ.txt", format = "fasta")
seqinr::consensus(seq_CG) 
nuc <- data.frame(x = (seqinr::consensus(seq_CG)))
as.matrix(nuc)
nuc$x
?seq.default
#slide function scrolls through variable set looking for matching pairs then outputs T/F if pair for the phrame is found.
CpG_Finder <- function(D, window, deslength){
  total <- length(D)
  lan <- deslength
  data <- D$x
  aa=lapply(x, function(data){
    x <- gregexpr("gc", data, perl = TRUE)
    res1 <- sum(attr(x[[1]], "match.length"))
    res1
  })
  bb=lapply(x, function(data){
    x <- gregexpr("g", data, perl = TRUE)
    res2 <- regmatches(data, x)
    res2
  })
  cc=lapply(x, function(data){
    x <- gregexpr("c", data, perl = TRUE)
    res3 <- regmatches(data, x)
    res3
  })
  
  result <- vector(length = length(data))
  for(i in 1:length(data)){result[i]}
  
  
  return(result)
}
CpG_Finder(nuc,2,500)


#### notes and comments section ####
# CpG graphing function()
n <- data.frame(read.csv("OverviewSelCoeff_BachelerFilter.csv"))

LvsF_CpG_Printer <- function(data_frame){
  maxnuc=1000
  plot(n$makesCpG[40:maxnuc],n$MeanFreq[40:maxnuc],
       log="y", ylab="frequency CpG non-CpG",cex.lab=1.3,
       xlab="Location",
       col="gray",t="n",pch=".", ylim=c(3.4*10^-4,1),xlim=c(40,979))
  
  points(c(n$MeanFreq[which(n$makesCpG=="1")],
           n$MeanFreq[which(n$makesCpG=="0")]),pch=p,col=co,
         bg=rgb(red=col2rgb(c)[1]/255,
                green=col2rgb(c)[2]/255,
                blue=col2rgb(c)[3]/255,
                maxColorValue = 1,alpha=0.8),
         cex=2)}#close for:
legpos=293;legposV=0.4
rect(legpos*3, 0.4*legposV, (legpos+42.5)*3, 1.7*legposV, density = NULL, angle = 45,col=alpha("white",1))


return()
}#close function

LvsF_CpG_Printer(n) #run function

#10-17-2017 nathan 12:53
#Formatted the group master script
# added sections