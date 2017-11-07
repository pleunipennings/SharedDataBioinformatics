#function to get WTAA

getWTAA<-function(DF){ #we make a fundtion that takes a data frame as input
    cons =  DF$wtnt #we are looking at the wtnt colum
    WTAA<-c() #make an empty vector WTAA to store the WT amino acids temporarily
    for(x in seq(1, length(cons) - 2, 3)){ # x is the looping variable, 
        #it goes from 1 until 2 less than the length of cons, and in 
        #goes in steps of three, 
        #so its value will be 1, 4, 7, 10 etc. 
        codon <- c(cons[x], cons[x+1], cons[x+2]) #get the codon that "belongs" 
        #to x, so for example, when x is 1, the codon is the nucleotides at 1, 2, 3
        new_AA <- seqinr::translate(codon) #translate the codon and store in new_AA
        WTAA[x] <- new_AA #put the Amino Acid in WTAA for x
        WTAA[x+1] <- new_AA #put the Amino Acid in WTAA for x+1
        WTAA[x+2] <- new_AA #put the Amino Acid in WTAA for x+2
    }
    if (length(which(names(DF)=="WTAA"))==0){ #check whether the dataframe 
        #has a WTAA column
        DF$WTAA=0} #If the WTAA column is not there, create it. 
    DF$WTAA<-WTAA #store the AA in the WTAA column
}