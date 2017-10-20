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
    
    

#### Mergerger ####

    ahub <- AnnotationHub()
    ahub["AH5086"]
    