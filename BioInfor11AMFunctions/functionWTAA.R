#function to get WTAA

getWTAA<-function(DF){
    cons =  DF$wtnt
    #cons<-c("A","A","C","C","G","T","C","G","T")
    WTAA<-c()
    for(x in seq(1, length(cons) - 2, 3)){
        codon <- c(cons[x], cons[x+1], cons[x+2])
        new_AA <- seqinr::translate(codon)
        WTAA[x] <- new_AA
        WTAA[x+1] <- new_AA
        WTAA[x+2] <- new_AA
    }
    if (length(which(names(DF)=="WTAA"))==0){
        DF$WTAA=0}
    DF$WTAA<-WTAA
}