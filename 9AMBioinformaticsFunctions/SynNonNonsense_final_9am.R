#9am
#function for syn/non/nonsence
#function for syn/non/nonsence
for (h in nrow(df)){
  if(df[h,"MUTAA"]== df[h,"WTAA"]){
    df[h,"TypeOfSite"] = "syn"
  }
  if(df[h,"MUTAA"] != df[h,"WTAA"]){
    if(df[h,"MUTAA"]=="*"){
      df[h,"TypeOfSite"] = "nonsense"}
    else {
      df[h,"TypeOfSite"] = "nonsyn"
    }
  }
}

