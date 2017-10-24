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
    
    

#### Mergerger ####
    
    HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
    consensus(HPIV1a)
    Error in obj[[1]]$tip.label : $ operator is invalid for atomic vectors
    ??read.alignment
    ??consensus
    example(ape::consensus)
    seqinr::consensus(HPIV1a)
    [1] "-" "g" "g" "a" "c" "a" "a" "g" "t" "c" "a" "c" "a" "g" "a" "c" "a" "t" "t"
    [20] "t" "g" "a" "t" "c" "t" "t" "a" "g" "t" "t" "a" "a" "a" "a" "c" "c" "t" "t"
    [39] "t" "a" "t" "a" "a" "t" "g" "g" "c" "t" "g" "g" "g" "c" "t" "a" "c" "t" "a"
    [58] "a" "g" "t" "a" "c" "t" "t" "t" "t" "g" "a" "c" "a" "c" "a" "t" "t" "c" "a"
    [77] "g" "c" "t" "c" "c" "a" "g" "g" "a" "g" "g" "a" "g" "t" "g" "a" "g" "a" "g"
    [96] "c" "a" "t" "c" "a" "a" "t" "a" "a" "g" "t" "c" "t" "g" "g" "c" "g" "g" "a"
    [115] "g" "g" "a" "g" "c" "a" "a" "t" "t" "a" "t" "a" "c" "c" "t" "g" "g" "t" "c"
    [134] "a" "a" "a" "g" "a" "a" "g" "t" "a" "c" "c" "g" "t" "t" "t" "c" "t" "g" "t"
    [153] "c" "t" "t" "c" "a" "c" "a" "t" "t" "a" "g" "g" "c" "c" "c" "g" "a" "g" "t"
    [172] "g" "t" "g" "a" "c" "a" "g" "a" "t" "g" "a" "t" "g" "c" "a" "g" "a" "t" "a"
    [191] "a" "a" "t" "t" "a" "t" "t" "a" "a" "t" "a" "g" "c" "a" "a" "c" "c" "a" "c"
    [210] "t" "t" "t" "c" "t" "t" "a" "g" "c" "c" "c" "a" "c" "t" "c" "a" "c" "t" "g"
    [229] "g" "a" "t" "a" "c" "a" "g" "a" "t" "a" "a" "a" "c" "a" "a" "c" "a" "c" "t"
    [248] "c" "t" "c" "a" "a" "a" "g" "a" "g" "g" "a" "g" "g" "a" "t" "t" "t" "t" "t"
    [267] "a" "g" "t" "a" "t" "c" "a" "c" "t" "c" "c" "t" "t" "g" "c" "a" "a" "t" "g"
    [286] "g" "c" "t" "t" "a" "c" "a" "g" "t" "a" "g" "c" "c" "c" "g" "g" "a" "g" "t"
    [305] "t" "a" "t" "a" "t" "c" "t" "c" "a" "c" "t" "a" "c" "a" "a" "a" "c" "g" "g"
    [324] "t" "g" "t" "c" "a" "a" "t" "g" "c" "t" "g" "a" "t" "g" "t" "c" "a" "a" "g"
    [343] "t" "a" "t" "g" "t" "g" "a" "t" "a" "t" "a" "c" "a" "g" "t" "a" "t" "a" "g"
    [362] "a" "g" "a" "g" "a" "g" "a" "t" "c" "c" "t" "a" "a" "a" "a" "g" "g" "a" "c"
    [381] "a" "a" "a" "a" "a" "c" "a" "g" "a" "t" "g" "g" "g" "t" "t" "c" "a" "t" "t"
    [400] "g" "t" "c" "a" "a" "a" "a" "c" "a" "a" "g" "a" "g" "a" "c" "a" "t" "g" "g"
    [419] "a" "g" "t" "a" "t" "g" "a" "a" "a" "g" "a" "a" "c" "a" "a" "c" "a" "g" "a"
    [438] "g" "t" "g" "g" "t" "t" "g" "t" "t" "c" "g" "g" "a" "c" "c" "t" "a" "t" "g"
    [457] "g" "t" "c" "a" "a" "c" "a" "a" "g" "a" "a" "c" "c" "c" "a" "t" "t" "g" "t"
    [476] "t" "c" "c" "a" "a" "g" "g" "g" "c" "a" "a" "a" "g" "a" "g" "a" "g" "a" "a"
    [495] "t" "g" "c" "g" "g" "a" "t" "c" "t" "a" "g" "a" "a" "g" "c" "a" "t" "t" "g"
    [514] "c" "t" "t" "c" "a" "g" "a" "c" "a" "t" "a" "t" "g" "g" "a" "t" "a" "t" "c"
    [533] "c" "t" "g" "c" "a" "t" "g" "t" "c" "t" "c" "g" "g" "a" "g" "c" "t" "a" "t"
    [552] "a" "a" "t" "a" "g" "t" "t" "c" "a" "a" "g" "t" "t" "t" "g" "g" "a" "t" "a"
    [571] "g" "t" "g" "c" "t" "g" "g" "t" "t" "a" "a" "g" "g" "c" "c" "a" "t" "a" "a"
    [590] "c" "a" "a" "g" "t" "a" "g" "t" "g" "c" "t" "g" "g" "t" "c" "t" "a" "a" "g"
    [609] "a" "a" "a" "a" "g" "g" "a" "t" "t" "c" "t" "t" "c" "a" "a" "t" "a" "g" "a"
    [628] "t" "t" "a" "g" "a" "a" "g" "c" "a" "t" "t" "c" "a" "g" "a" "c" "a" "a" "g"
    [647] "a" "t" "g" "g" "a" "a" "c" "c" "g" "t" "c" "a" "a" "a" "a" "g" "t" "g" "c"
    [666] "t" "c" "t" "g" "g" "t" "c" "t" "t" "c" "a" "c" "a" "g" "g" "a" "g" "a" "c"
    [685] "a" "c" "a" "g" "t" "t" "g" "a" "a" "g" "g" "c" "a" "t" "t" "g" "g" "t" "g"
    [704] "c" "a" "g" "t" "g" "a" "t" "g" "a" "g" "g" "t" "c" "a" "c" "a" "a" "c" "a"
    [723] "a" "a" "g" "c" "t" "t" "a" "g" "t" "a" "t" "c" "t" "c" "t" "t" "a" "t" "g"
    [742] "g" "t" "a" "g" "a" "g" "a" "c" "t" "c" "t" "a" "g" "t" "g" "a" "c" "t" "a"
    [761] "t" "g" "a" "a" "c" "a" "c" "a" "t" "c" "c" "a" "g" "g" "t" "c" "a" "g" "a"
    [780] "t" "c" "t" "a" "a" "c" "t" "a" "c" "a" "t" "t" "a" "g" "a" "g" "a" "a" "g"
    [799] "a" "a" "c" "a" "t" "t" "c" "a" "g" "a" "t" "t" "g" "t" "a" "g" "g" "a" "a"
    [818] "a" "t" "t" "a" "c" "a" "t" "a" "a" "g" "a" "g" "a" "t" "g" "c" "a" "g" "g"
    [837] "a" "t" "t" "a" "g" "c" "a" "t" "c" "t" "t" "t" "c" "a" "t" "g" "a" "a" "c"
    [856] "a" "c" "c" "a" "t" "c" "a" "a" "g" "t" "a" "t" "g" "g" "t" "g" "t" "a" "g"
    [875] "a" "a" "a" "c" "g" "a" "a" "g" "a" "t" "g" "g" "c" "c" "g" "c" "c" "t" "t"
    [894] "g" "a" "c" "a" "c" "t" "a" "t" "c" "a" "a" "a" "c" "c" "t" "g" "a" "g" "a"
    [913] "c" "c" "a" "g" "a" "t" "c" "t" "a" "a" "a" "c" "a" "a" "a" "c" "t" "g" "a"
    [932] "g" "a" "a" "g" "c" "c" "t" "t" "g" "t" "t" "g" "a" "t" "a" "t" "c" "t" "a"
    [951] "t" "c" "t" "a" "t" "c" "a" "a" "a" "g" "g" "g" "a" "g" "c" "c" "c" "g" "a"
    [970] "g" "c" "c" "c" "c" "t" "t" "t" "t" "a" "t" "a" "t" "g" "t" "a" "t" "a" "c"
    [989] "t" "c" "a" "g" "a" "g" "a" "c" "c" "c" "a" "g"
    [ reached getOption("max.print") -- omitted 14473 entries ]
    > length(seqinr::consensus(HPIV1a))
    [1] 15473
    
  #  1/ make consensus data of our sample DNA (for one file only!)
n <- data.frame(Pos = c(1:length(seqinr::consensus(HPIV1a))),
                WTnt = (seqinr::consensus(HPIV1a)),
                Trans = c("frequency()")
                )

alphabetFrequency(n$WTnt, c("c","g","a","t"), as.prob = T, baseOnly=TRUE
                  )


consensusString(HPIV1a, baseOnly=TRUE, ambiguityMap=IUPAC_CODE_MAP)
head(n)
tail(n)
  #  2/ determine the frequency of mutation at each base
  #  3/ create a function that will plot the cpg /posible cpg





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
################### new section of shame ###################