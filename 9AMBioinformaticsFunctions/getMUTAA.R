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
