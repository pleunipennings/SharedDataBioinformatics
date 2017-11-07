
library("Biostrings")
fasta2dataframe=function(seqdump.txt)
    {
    s = readDNAStringSet(seqdump.txt)
    RefSeqID = names(s)
    RefSeqID = sub(" .*", "", RefSeqID) 
   
    #erase all characters after the first space: regular expression matches a space followed by any sequence of characters and sub replaces that with a string having zero  characters 
    
    for (i in 1:length(s))
        
        {
        seq[i]=toString(s[i])
    }
    
    RefSeqID_seq=data.frame(RefSeqID,seq)
    return(RefSeqID_seq)
}
RefSeqID



