#  Note that this depends on the columns MUTAA and WTAA being already complete!

functionSynNonSyn<-function(DF){
<<<<<<< HEAD
    if (length(which(names(DF))=="MUTAA")==0){
        print("Oh oh there is a problem. No MUTAA column!")
        return(0)}
=======
    if (length(which(names(DF)=="MUTAA"))==0){
        print("Oh oh there is a problem. No MUTAA column!")
        return(0)}
    if (length(which(names(DF)=="TypeOfSite"))==0){
        DF$TypeOfSite=0}
>>>>>>> 92263dd95301c435bc46a2d1fe1055acc8d2ce69
    for (h in 1:nrow(DF)){
        if(DF$MUTAA[h]== DF$WTAA[h]){
            #   if(DF[h,"MUTAA"]== DF[h,"WTAA"]){
            DF[h,"TypeOfSite"] = "syn"
        }
        if(DF[h,"MUTAA"] != DF[h,"WTAA"]){
            if(DF[h,"MUTAA"]=="*"){
                DF[h,"TypeOfSite"] = "nonsense"}
            else {
                DF[h,"TypeOfSite"] = "nonsyn"
            }
        }
    }
    DF$TypeOfSite<-as.factor(DF$TypeOfSite)
}
<<<<<<< HEAD
=======

>>>>>>> 92263dd95301c435bc46a2d1fe1055acc8d2ce69
