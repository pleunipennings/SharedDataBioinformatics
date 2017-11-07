#  Note that this depends on the columns MUTAA and WTAA being already complete!

<<<<<<< HEAD
functionSynNonSyn<-function(dengue_df){
    if (length(which(names(dengue_df))=="MUTAA")==0){
        print("Oh oh there is a problem. No MUTAA column!")
        return(0)}
    for (h in 1:nrow(dengue_df)){
        if(dengue_df$MUTAA[h]== dengue_df$WTAA[h]){
            #   if(dengue_df[h,"MUTAA"]== dengue_df[h,"WTAA"]){
            dengue_df[h,"TypeOfSite"] = "syn"
=======
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
>>>>>>> d0aea278e91ebd92b73d53de3c019f40c51ad8aa
        }
        if(dengue_df[h,"MUTAA"] != dengue_df[h,"WTAA"]){
            if(dengue_df[h,"MUTAA"]=="*"){
                dengue_df[h,"TypeOfSite"] = "nonsense"}
            else {
                dengue_df[h,"TypeOfSite"] = "nonsyn"
            }
        }
    }
    dengue_df$TypeOfSite<-as.factor(dengue_df$TypeOfSite)
}
<<<<<<< HEAD
=======

>>>>>>> 92263dd95301c435bc46a2d1fe1055acc8d2ce69
