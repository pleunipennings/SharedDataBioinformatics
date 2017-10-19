#### Master Script ####

#TEST 

#Fig 2: location vs frequency CpG non-CpG.
  #Plot num column vs MeanFreq column in a scatterplot. The points should be colored
  #depending on whether they are Cpg / not CpG.

#remember [PULL > Comment > Commit > PUSH]




#### function to find CpG Islands ####

CpG_finder <- function(data, window, step){
        total <- length(data)
        spots <- seq(from = 1, to = (total-window), by = step)
        result <- vector(length = length(spots))
        for(i in 1:length(spots)){result[i]}
        return(result)
    }
slide_function(seqT,3,500)


####NEW CODE FROM GITHUB

CpG-Islander
-============
    
    -CpG Islander
+#
    +# This script is designed for makeCGI package installation and execution.
    +#
    +# The papers about CGI search using HMM can be found here:
    +# http://www.ncbi.nlm.nih.gov/pubmed/20212320
    +# http://www.ncbi.nlm.nih.gov/pubmed/19777308
    +#
    +#
    +# Installation process
    +# 1. Install "BSgenome" and "Biostrings" Bioconductor packages (if you are in trouble see the instructions below)
    +
    +# Istalling Bioconductor packages
    +#
    +# IF YOUR INTERNET CONNECTION USES PROXY
    +# Replace username, password, proxy ip and port and run following in the terminal:
    +# export http_proxy=http://username:password@proxy_address:proxy_port/
    +# export HTTP_PROXY=http://username:password@proxy_address:proxy_port/
    +#
    +# If you are using RStudio then restart it
    +#
    +# Set proxy and check it
    +#Sys.setenv(http_proxy="http://username:password@proxy_address:proxy_port/")
    +# Check that proxy set up succeed
    +#Sys.getenv("http_proxy")
    +
    +# INSTALL BIOCONDUCTOR PACKAGES
    +source("http://bioconductor.org/biocLite.R")
+biocLite()
+#
    +# Install makeCGI dependencies
    +biocLite("BSgenome")
+biocLite("Biostrings")
+
    +# 2. Download and install makeCGI from the following page:
    +# http://rafalab.jhsph.edu/CGI/
    +# You can install it via Tools/Install packages menu of RStudio or run the following commands
    +makecgi_directory="changeme"
+setwd(makecgi_directory)
+install.packages("makeCGI_1.2.tar.gz", repos = NULL, type = "source")

#### function to plot CpG Islands ####
CpG_plotter <- function(){
    
}

#### notes and comments section ####
#10-17-2017 nathan 12:53
    #Formatted the group master script
    # added sections
# pure Fasta format Upload
HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
HPIV1b = read.alignment("humanparainfluenzavirus1_F.fasta_pruned.mu.trim05", format = "fasta")
HPIV1c = read.alignment("humanparainfluenzavirus1_HN.fasta_pruned.mu.trim05", format = "fasta")
HPIV3a = read.alignment("humanparainfluenzavirus3.fasta_pruned.mu.trim05", format = "fasta")
HPIV3b = read.alignment("humanparainfluenzavirus3_HN.fasta_pruned.mu.trim05", format = "fasta")
#average size of each sample is 15,500~ approx
# none of this worked!
#       str(HPIV1a)
#       as.table(HPIV1a, stringsAsFactors=FALSE)
#       M <- as.table(read.alignment("clean_HPIV1.txt", format = "fasta"))
#       M

#       B <-read.csv("Clean_HPIV1.csv")
#       B
#       str(B)

# I took each fasta file and converteded it into a csv tabulated format
# the T stands for Tabulated format or .csv
HPIV1aT = read.csv("HPIV1aT.csv") #from HPIV1a = read.alignment("humanparainfluenzavirus1.fasta_pruned.mu.trim05", format = "fasta")
HPIV1bT = read.csv("HPIV1bT.csv") #from HPIV1b = read.alignment("humanparainfluenzavirus1_F.fasta_pruned.mu.trim05", format = "fasta")
HPIV1cT = read.csv("HPIV1cT.csv") #from HPIV1c = read.alignment("humanparainfluenzavirus1_HN.fasta_pruned.mu.trim05", format = "fasta")
HPIV3aT = read.csv("HPIV3aT.csv") #from HPIV3a = read.alignment("humanparainfluenzavirus3.fasta_pruned.mu.trim05", format = "fasta")
HPIV3bT = read.csv("HPIV3bT.csv") #from HPIV3b = read.alignment("humanparainfluenzavirus3_HN.fasta_pruned.mu.trim05", format = "fasta")

#average size of each sample is 15,500~ approx
str(HPIV1aT)
head(HPIV1aT)





