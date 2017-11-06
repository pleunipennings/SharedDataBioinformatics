#following packages are required >>
library(graphics)
library(seqinr)

####### CpG graphing function() #########

LvsF_CpG_Printer <- function(n){ # value "n" will represent our data.frame of use
    if (T) {n$makesCpG <- n$makesCpG+1} #this adds one value to the "0" and "1" boolian values to "1" and "2" this is for the as.integer condition to work
    YCpG <- which(n$makesCpG=="2") #lists which variables return a "2" these make a CpG island when mutated "yes cpg or Y"
    NCpG <- which(n$makesCpG=="1") #lists which variables return a "1" these do not make a CpG island when mutated "No cpg or N"
    x1 <- n$MeanFreq[YCpG] #create a Value that looks at mean frequency by yes cpg
    x2 <- n$MeanFreq[NCpG] #create a Value that looks at mean frequency by no cpg
    plot.default(x = c(x1, x2), 
                 xlab = "Location", ylab = "Frequency", main = "Location vs frequency CpG non-CpG Graph",
                 col = (as.integer(n$makesCpG)),
                 log = "y"
    )#close plot.default
    return()#close return
}#close function

LvsF_CpG_Printer(n) #run function

