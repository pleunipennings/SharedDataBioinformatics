#Team: Livia, Gabriella, Deshawn


###Big AA change###
big_aa_change <- function(df){
    
    #Function for amino acid categories
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
    
    #Create empty vectors
    WTAAcategory <- c()
    MUTAAcategory <- c()
    bigAAchange <- c()
    WTAA<-df$WTAA
    MUTAA<-df$MUTAA
    
    #Assign wild type AA category
    for(j in 1:nrow(df)){
        WTAAcategory[j]=amCat(WTAA[j])
    }
    
    #Assign mutated AA category
    for(j in 1:nrow(df)){
        MUTAAcategory[j]=amCat(MUTAA[j])
    }
    
    #Loop for drastic change or not 
    for(i in 1:nrow(df)){
        if (WTAAcategory[i]==MUTAAcategory[i]){
            bigAAchange[i]= "0"
        }
        if (WTAAcategory[i]!=MUTAAcategory[i]){
            bigAAchange[i] = "1"
        }
    }
    
    #Create "bigAAchange" column if not already
    if (length(which(names(df)=="bigAAchange"))==0){
        df$bigAAchange=0}
    
    #Inserting value into column
    df$bigAAchange <- bigAAchange
    
    #return data frame
    return(df)
}