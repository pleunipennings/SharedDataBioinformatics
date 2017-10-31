#  Note that this depends on the columns MUTAA and WTAA being already complete!

for (h in nrow(df)){
  
  if(df[h,"MUTAA"]== df[h,"WTAA"]){
    
    df[h,"TypeOfSite"] = "syn"
    
  }
  
  if(df[h,"MUTAA"] != df[h,"WTAA"]){
    
    if(df[h,"MUTAA"]=="*"){
      
      df[h,"TypeOfSite"] = "non-sense"}
    
    else {
      
      df[h,"TypeOfSite"] = "non-syn"
      
    }
    
  }
  
}