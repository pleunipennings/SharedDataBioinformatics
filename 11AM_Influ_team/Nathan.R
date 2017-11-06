#Nathan
if (TRUE) {
    library(ggplot2)
    library(plyr)
    library(grid)
    library(scales)
    library(gridExtra)
    library(seqinr)
}

#### group function for submission ####
HPIV1a
seqinr::consensus(HPIV1a)
length(seqinr::consensus(HPIV1a))
DF <- data.frame(WTAA = (seqinr::consensus(HPIV1a)), 
                 MUTAA(), 
                 TypeOfSite(),
                 MakesCpG(), 
                 bigAAchange()
)

HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
Pos<-c(1:15473)
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

head(HPIVdf)
tail(HPIVdf)
HPIVdf$MeanFreq

HPIVdf$WtNt<- seqinr::consensus(HPIV1a)
#current function as is 
# need to get legend to work
####### CpG graphing function() #########

LvsF_CpG_Printer <- function(n){ # value "n" will represent our data.frame of use
    if (T) {n$makesCpG <- n$makesCpG+1} #this adds one value in "n" data.frame to the "0" and "1" boolian values to "1" and "2" this is for the as.integer condition to work
    YCpG <- which(n$makesCpG=="2") #lists which variables return a "2" these make a CpG island when mutated "yes cpg or Y"
    NCpG <- which(n$makesCpG=="1") #lists which variables return a "1" these do not make a CpG island when mutated "No cpg or N"
    x1 <- n$MeanFreq[YCpG] #create a Value that looks at mean frequency by yes cpg
    x2 <- n$MeanFreq[NCpG] #create a Value that looks at mean frequency by no cpg
    plot.default(x = c(x1, x2), #y values not listed this is intended
                 xlab = "Position Location", ylab = "Frequency", main = "Location vs frequency CpG non-CpG Graph",
                 col = (as.integer(n$makesCpG)), #colors the values by their integer, using default association. 
                 log = "y" #the plot is put into a log scale 
    )#close plot.default
    legend('topright', #add legend in top right corner
           legend = c(unique(n$makesCpG)), #names of each category based on the factors in n$makesCpG
           col=unique(as.integer(n$makesCpG)), #colors the factors as is found by integers in n$makeCpG
           pch=21 #symbols matching dataframe's factors
           ) #close legend
  
    return()# prints out plot with legend 
}#close function

LvsF_CpG_Printer(n) #run function

NULL
?month.abb
legend("topright", 
       inset= 0.01, #inset legend off from border
       c("YesCpG", "NoCpG"),
       legend = c(makesCpG == 2, makesCpG == 1), #names of each category based on factors, alphabetical order of category 1-5 of TypeOfSites
       pch= (1), 
       col= c(as.integer(n$makesCpG)),  #colors matching dataframe's factors 1:5 of DF$TypeOfSite
       cex=.75,   #specify scale of whole box of legend to not block data
       pt.cex=3,  #specify point size in legend
       #remove legend border if "n" is specificed. "o" displays border
       bty="o",
       #specific box border thickness/width
       box.lwd=2,
       #text.width change
       text.width=120,
       #justify text legend, xjust=0 is left justified, xjust=1 means right justified
       xjust=0
)

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

functionSynNonSyn <- function(DF){
    if(length(which(names(DF))=="MUTAA")==0)
for (h in 1:nrow(DF)){
    if(DF$MUTAA[h]== DF$WTAA[h]){
   #if(DF[h,"MUTAA"]== DF[h,"WTAA"]){
        DF[h, "TypeOfSite"] = "syn"
    }
    if(df)
    }
}




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
###### CPG finder ######
CpG_finder <- function(new_virus_data){
    #reads data into function as CSV file
    virus_data <- read.csv(new_virus_data)
    
    #singles out column of nucleotides -- this assumes that the new datafile has the same column headers such as WTnt as the HIV file does
    WTnt <- virus_data$WTnt
    
    #gets the length of the data file for use in a loop
    data_length <- nrow(virus_data)
    
    #creates an empty vector (like a list) of the same length as the data file to be used to record the results of the loop below
    makesCpG <- vector(mode = "numeric", length = data_length)
    
    #loop that determines if a CpG could occur due to mutation at each spot in the list of nucleotides (WTnt)
    #loops from row 1 to the the last row of the column WTnt
    for(x in 1:data_length){
        
        #assigns a name (current_nuc) to the nucleotide at row x in WTnt and makes the nucleotide capitalized, in case the data uses lower case letters
        current_nuc <- toupper(WTnt[x])
        
        #assigns a name (current_neighbor) to the nucleotide in the next row down in WTnt and makes the nucleotide capitalized, in case the data uses lower case letters
        current_neighbor <- toupper(WTnt[x + 1])
        
        #begins a loop that compares the two nucleotides that were just isolated
        #this if statement is only evaluated if the current nucleotide is a T
        if(current_nuc == "T"){
            #this if statement is only evaluated if the neighbor is a G
            if(current_neighbor == "G"){
                #if the current nucleotide is a T and the neighbor is a G (thus a mutation can cause a CpG), then the of the current nucleotide in the list of 0s we made is changed to a 1
                makesCpG[x] = 1
            }
            #otherwise, nothing is changed
            else{}
        }
        #repeat the same loop as before, but modified slightly to look for CA pairs
        if(current_nuc == "C"){
            if(current_neighbor == "A"){
                #here, instead of changing the position of the current nuc to a 1, we change the neighbor's position because that is where the mutation would be, hence x+1
                makesCpG[x+1] = 1
            }
            else{}
        }
        else{}
        #print(c("Data length is ", nrow(virus_data), "List length is ", length(makesCpG)))
    }
    
    #append the list of 0s and 1s we created to the end of the data file we imported
    virus_data$makesCpG <- makesCpG
    
    #when the entire function is done, return the amended data file
    return(virus_data)
}

##### Function that makes synonymous/nonsynonymous column for data frame ######
#  Note that this depends on the columns MUTAA and WTAA being already complete!

functionSynNonSyn<-function(DF){
    if (length(which(names(DF))=="MUTAA")==0){
        print("Oh oh there is a problem. No MUTAA column!")
        return(0)}
    for (h in 1:nrow(DF)){
        if(DF$MUTAA[h]== DF$WTAA[h]){
            #   if(DF[h,"MUTAA"]== DF[h,"WTAA"]){
            DF[h,"TypeOfSite"] = "syn"
        }
        if(DF[h,"MUTAA"] != DF[h,"WTAA"]){
            if(DF[h,"MUTAA"]=="*"){
                DF[h,"TypeOfSite"] = "nonsense"}
            else {
                DF[h,"TypeOfSite"] = "nonsyn"
            }
        }
    }
    DF$TypeOfSite<-as.factor(DF$TypeOfSite)
}

##### Big AA CHANGE  #####

#Read in shortened entero data
HPIVseqs<-read.fasta("humanparainfluenzavirus1.fasta_pruned.mu.trim05")

#Align entero data
HPIValigned<-read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format="fasta")

#Gets DNA into matrix form
HPIVmatrix<-read.dna("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format="fasta", as.character = TRUE)

#Gets WTnt for each nt position
cons<-seqinr::consensus(HPIValigned)

#Insert concensus into dataframe
HPIVdf$WtNt=cons
WTAA<- NULL
MUTAA<-NULL
#removes list data

#Function to return transition mutation
transition<-function(basepair){
    #basepair<-("A", "C", "T", "G"),
    if(basepair=="a") {return("g")}
    if(basepair=="g") {return ("a")}
    if(basepair=="t") {return ("c")}
    if(basepair=="c") {return ("t")}
}
head(HPIVdf)
#Loop to insert the freq of transition mutations
for(i in 1:nrow(HPIVdf)){
    HPIVdf$MeanFreq[i]<-length(which(HPIVmatrix[,i]==transition(cons[i])))/
        (nrow(HPIVmatrix))
}

#for loop for transitioning each position and translating it 
for(x in seq(1, length(cons), 3)){
    codon <- c(cons[x], cons[x+1], cons[x+2])
    mutated_codon <- codon
    if(codon[1] == "a"){
        mutated_codon <- replace(x=mutated_codon, values=c("g", codon[2], codon[3]))
    }
    if(codon[1] == "g"){
        mutated_codon <- replace(x=mutated_codon, values=c("a", codon[2], codon[3]))
    }
    if(codon[1] == "c"){
        mutated_codon <- replace(x=mutated_codon, values=c("t", codon[2], codon[3]))
    }
    if(codon[1] == "t"){
        mutated_codon <- replace(x=mutated_codon, values=c("c", codon[2], codon[3]))
    }
    MUTAA[x] <- translate(mutated_codon)
    
    if(codon[2] == "a"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "g", codon[3]))
    }
    if(codon[2] == "g"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "a", codon[3]))
    }
    if(codon[2] == "c"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "t", codon[3]))
    }
    if(codon[2] == "t"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], "c", codon[3]))
    }
    MUTAA[x+1] <- translate(mutated_codon)
    
    if(codon[3] == "a"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "g"))
    }
    if(codon[3] == "g"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "a"))
    }
    if(codon[3] == "c"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "t"))
    }
    if(codon[3] == "t"){
        mutated_codon <- replace(x=mutated_codon, values=c(codon[1], codon[2], "c"))
    }
    MUTAA[x+2] <- translate(mutated_codon)
}

for(x in seq(1, length(cons) - 2, 3)){
    codon <- c(cons[x], cons[x+1], cons[x+2])
    new_AA <- translate(codon)
    WTAA[x] <- new_AA
    WTAA[x+1] <- new_AA
    WTAA[x+2] <- new_AA
}

#places respective things in columns
HPIVdf$WTAA<- WTAA
HPIVdf$MUTAA <- MUTAA

#Amino Acid Changes 
pos <- "R|H|K"
neg <- "D|E"
unc <- "S|T|N|Q"
spe <- "C|U|G|P"
hyd <- "A|I|L|F|M|W|Y|V"
amCat <- function(AA){
    if(regexpr(pos, AA) > 0){ return(0) }
    if(regexpr(neg, AA) > 0){ return(1) }
    if(regexpr(unc, AA) > 0){ return(2) }
    if(regexpr(spe, AA) > 0){ return(3) }
    if(regexpr(hyd, AA) > 0){ return(4) }
    return(5)
}

#Assign wild type AA category
for(j in 1:nrow(HPIVdf)){
    HPIVdf$WTAAcat[j]=amCat(HPIVdf$WTAA[j])
}

#Assign mutated AA category
for(j in 1:nrow(HPIVdf)){
    HPIVdf$MUTAAcat[j]=amCat(HPIVdf$MUTAA[j])
}

#Loop for drastic change or not 
for(i in 1:nrow(HPIVdf)){
    if (HPIVdf$WTAAcat[i]==HPIVdf$MUTAAcat[i]){
        HPIVdf$bigAAchange[i]= "0"
    }
    if (HPIVdf$WTAAcat[i]!=HPIVdf$MUTAAcat[i]){
        HPIVdf$bigAAchange[i] = "1"
    }
}

######### Code for Figure 2 ###############

summary(HPIVdf$MeanFreq)
summary(log(HPIVdf$MeanFreq))
plot(
    #x vector
    HPIVdf$Pos,
    #y vector
    HPIVdf$TrNtFreq,
    #make black empty circles as symbol
    pch=21,
    #make outline of symbol black
    col= "black",
    #fill inside of point with color by factor category "TypeOfSite" bg=
    bg=HPIVdf$TypeOfSite,
    #Title label
    main = "HIV Practice Data",
    #x axis label
    xlab ="Location on Sequence", 
    #y axis label
    ylab ="Mean Frequency of Mutation",
    #cex changes point size
    cex=2,
    #grid superimposes grid onto plot 
    #nx and ny describes x and y axis 
    #NA will automatically set x-axis to default plot ticks 
    grid(nx = NA, ny = NULL, col = "black", lty = "dotted",
         lwd = par("lwd"), equilogs = TRUE),
    #supress y axis drawing by plot fxn, put # in front to not supress
    yaxt="n",
    # or do log of y axis, delete # symbol
    log="y"
)

#axis function to write y axis with only scale by 10s. 
#USE YOUR OWN APPROPRIATE NUMBERS
axis(
    #2 is to specify left axis (aka y axis)
    2,
    at=c(0.0001,0.001,0.01,0.1)
    # dont need to define labels if same as numbers on axis
    #labels=aty
    #tck marks
    #tck=-0.01
)

#add legend in top right corner
legend("topright", 
       #inset legend off from border
       inset= 0.01,
       #names of each category based on factors, alphabetical order of category 1-5 of TypeOfSites
       legend = levels(HPIVdf$TypeOfSite), 
       #symbols matching dataframe's factors 1:5 of data$TypeOfSite
       pch=21,
       #colors matching dataframe's factors 1:5 of data$TypeOfSite
       col="black",
       #fill colors of circle matching points of plot
       pt.bg=c(1:5),
       #specify scale of whole box of legend to not block data
       cex=.75,
       #specify point size in legend
       pt.cex=3,
       #remove legend border if "n" is specificed. "o" displays border
       bty="o",
       #specific box border thickness/width
       box.lwd=2,
       #specify box border type
       #box.lty=,
       #text.width change
       text.width=10,
       #justify text legend, xjust=0 is left justified, xjust=1 means right justified
       xjust=1
)


######## code for fig 3 #########
### THESE PACKAGES MUST BE INSTALLED ONLY ONCE:
install.packages("ggplot2")
install.packages("plyr")
install.packages("grid")
install.packages("ape")
install.packages("seqinr")
install.packages("gridExtra")

### SET YOUR WORKING DIRECTORY BEFORE STARTING 
### AND READ CSV FILE (WHATEVER YOU CALLED YOUR PRACTICE DATA):

setwd("~/Desktop/Bioinformatics (Biol 738)/Midterm Project/")
read.csv("HIV practice data.csv") -> CSV

library(ggplot2)
library(plyr)
library(grid)
library(scales)
library(gridExtra)

#### FUNCTION FOR FIGURE 3 STARTS HERE

pug<-function(dfx){
    
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
    #all green points for the right NON-synonomous site graph
    g <- frequenciesOfNONSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "a" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
    k <- frequenciesOfNONSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "t" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
    o <- frequenciesOfNONSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "c" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
    q <- frequenciesOfNONSynAmutsNONCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "g" & (CSV$bigAAChange == "0") &(CSV$makesCpG == "0"))),"MeanFreq"]
    
    #blue dots of right graph that show a CpG-forming mutation
    h <- frequenciesOfNONSynAmutsCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "a") & (CSV$bigAAChange == "0") & (CSV$makesCpG == "1")),"MeanFreq"]
    l <- frequenciesOfNONSynTmutsCP <- CSV[which(((CSV$TypeOfSite == "nonsyn"  )) & (CSV$WTnt == "t") & (CSV$bigAAChange == "0") & (CSV$makesCpG == "1")),"MeanFreq"]
    
    #all yellow points for the right NON-SYN site graph
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
        
        geom_point( data = dfval, mapping = aes(x="a", y= val_a), colour = "green", size = 0.1) +xlab("Mutation Type") + ylab("Mutation Frequency") +
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
    
    
    
    
    
    
    # Creating Non-Synonymous Plot
    NonSynplt <- ggplot() +
        
        geom_point( data = dfval, mapping = aes(x="g", y= val_g), colour = "green", size = 0.1) +xlab("Mutation Type") + ylab("Mutation Frequency") +
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
                    Fig. 3: Frequencies", gp=gpar(fontface="bold", fontsize = 16, cex = 1))
    grid.arrange( top =title1,Synplt + ggtitle('Synonymous Sites'), NonSynplt + ggtitle('Non-synonymous Sites'),  nrow=1)
    
    # Creating the plot keys
    grid.text("KEY:
              Green = No Drastic AA change (non-CpG forming)
              Blue = No Drastic AA change (CpG forming)
              Orange = Drastic AA change (non-CpG-forming)
              Red = Drastic AA change (CpG forming) ", 
              x = unit(2, "cm"), y = unit(0.25,"cm"), just = "left", vjust = unit(0.0,"cm"))
    
    
    #sink("Routput(loop).txt")
    #print(dfval)
    #print(dfstats)
    #class(dfstats)
    #sink()
    
    
}

### ENTER IN WHATEVER YOU NAMED THE PRACTICE DATA ###
anything<-read.csv("HIV practice data.csv")
pug(anything)
#######

source("functionWTAA.R")
