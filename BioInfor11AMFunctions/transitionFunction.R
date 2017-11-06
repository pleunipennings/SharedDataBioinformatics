

#Function to return transition mutation
transition<-function(basepair){
    #basepair<-("A", "C", "T", "G"),
    if(basepair=="a") {return("g")}
    if(basepair=="g") {return ("a")}
    if(basepair=="t") {return ("c")}
    if(basepair=="c") {return ("t")}
    if(basepair=="-") {return ("-")}
}