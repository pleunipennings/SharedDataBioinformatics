#  Note that this depends on the columns MUTAA and WTAA being already complete!

functionSynNonSyn<-function(dengue_df){

    #Create "TypeOfSite" column if not already
    if (length(which(names(dengue_df)=="TypeOfSite"))==0){
        dengue_df$TypeOfSite<-0}
    
 
    
    
    if (length(which(names(dengue_df)=="MUTAA"))==0){
        print("Oh oh there is a problem. No MUTAA column!")
        return(0)}
    for (h in 1:nrow(dengue_df)){
        if(dengue_df$MUTAA[h]== dengue_df$WTAA[h]){
            #   if(dengue_df[h,"MUTAA"]== dengue_df[h,"WTAA"]){
            dengue_df[h,"TypeOfSite"] = "syn"}
        if(dengue_df[h,"MUTAA"] != dengue_df[h,"WTAA"]){
            if(dengue_df[h,"MUTAA"]=="*"){
                dengue_df[h,"TypeOfSite"] = "nonsense"}
            else {
                dengue_df[h,"TypeOfSite"] = "nonsyn"
            }
        }
    }
    dengue_df$TypeOfSite<-as.factor(dengue_df$TypeOfSite)
    return(dengue_df)
    }
