# This function determines if the mutated amino acid from a transition mutation is a
# big amino acid change or not.

#This function will first require the use of two other functions (all are below).

# the input is a datafame with a "wtnt" column



###Get Wildtype amino acid###
getWTAA<-function(df){
  
  #Assign consensus to a variable
  cons =  df$wtnt
  
  #Create empty vector for wildtype amino acid
  WTAA <- c()
  
  #Loop for translating consensus
  for(x in seq(1, length(cons) - 2, 3)){
    codon <- c(cons[x], cons[x+1], cons[x+2])
    new_AA <- seqinr::translate(codon)
    WTAA[x] <- new_AA
    WTAA[x+1] <- new_AA
    WTAA[x+2] <- new_AA
  }
  
  #Create "WTAA" column if not already
  if (length(which(names(df)=="WTAA"))==0){
    df$WTAA=0}
  
  #Insert value into column
  df$WTAA<-WTAA
  
  #Return the data frame
  return(df)
}



###Get mutated amino acid###
getMUTAA <- function(df){
  
  #Assign consensus to a variable
  cons = df$wtnt
  
  #Create empty vector for mutated amino acid
  MUTAA <- c()
  
  #Loop for mutated codon
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
  
  #Create "MUTAA" category if not already
  if (length(which(names(df)=="MUTAA"))==0){
    df$MUTAA=0}
  
  #Insert value into column
  df$MUTAA<-MUTAA
  
  #Return data frame
  return(df)
}



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

#Team: Livia, Gabriella, Deshawn