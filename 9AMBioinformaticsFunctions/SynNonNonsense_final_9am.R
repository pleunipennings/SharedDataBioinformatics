#9am
#function for syn/non/nonsence
#function for syn/non/nonsence
for (h in nrow(df)){ #looks at each row in the dataframe df
  if(df[h,"MUTAA"]== df[h,"WTAA"]){ #if the value in the MUTAA column for the row of interest is equal to the value in the WTAA of the same row
    df[h,"TypeOfSite"] = "syn" #then insert "syn" for the value in the TypeOfSite column for this row
  }
  if(df[h,"MUTAA"] != df[h,"WTAA"]){ #if the value in the MUTAA column for the row of interest is NOT equal to the value in the WTAA of the same row
    if(df[h,"MUTAA"]=="*"){ #and it the MUTAA value equal "*"
      df[h,"TypeOfSite"] = "nonsense"} #then insert "nonsense" for the value in the TypeOfSite column for this row
    else {
      df[h,"TypeOfSite"] = "nonsyn" #else insert "nonsyn" for the value in the TypeOfSite column for this row
    }
  }
}
#RM
