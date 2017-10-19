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
    spots <- seq(from = 1, to = (total-window), by = step)
    result <- vector(length = length(spots))
    for(i in 1:length(spots)){result[i]}
    return(result)
}
slide_function(seqT,3,500)

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
    