#### Plot CpG By Feq & Location ####

# written and coded by 11am Human parainfluenza virus team:
    # amir J.
    # annie S.
    # andrew M.
    # nathan O.

# outside contributors:
    # alburt

# description of function:
    # LvsF_CpG_Printer() functions as its name implies. 
    # the function when run looks at which values that make a cpg vs those that dont.
    # then it looks at the  mean freqency of CpG formation at these values in the dataframe.
    # this is then plotted in a log scale.
    # values that make CpG's are listed as 2
    # values that do not make CpG's are listed as 1
    # these are color coded as well!

##### following packages are required #####

# Ctrl run the following if statement
if (TRUE) {
    library(ggplot2)
    library(plyr)
    library(grid)
    library(scales)
#    library(gridExtra)
    library(seqinr)
}
##### assumptions #####
    # above packages are installed and loaded
    # data.frame has the value $MeanFreq
    # data.frame has the value $makesCpG
    # values in $makesCpG are listed as "0" and "1"
        # If not listed as 0's or 1's 
        #  or
        # If title is wrong prints out: 
          # Error in `$<-.data.frame`(`*tmp*`, "MakesCpG", value = numeric(0)) 

####### function() for plotting location by frequency for CpG's occurence  #########

LvsF_CpG_Printer <- function(n){ # value "n" will represent our data.frame of use

    #    if (T) {n$makesCpG <- n$makesCpG+1} #this adds one value in "n" data.frame to the "0" and "1" boolian values to "1" and "2" this is for the as.integer condition to work
    #PSP: instead, let's use n$makesCpG+1 directly for the colors

    smallestvalue = min(n$MeanFreq[n$MeanFreq>0]) # a small value to add to all frequencies to be able to see zeros on a log scale
    
    YCpG <- which(n$makesCpG==1) #lists which variables return a "2" these make a CpG island when mutated "yes cpg or Y"
    NCpG <- which(n$makesCpG==0) #lists which variables return a "1" these do not make a CpG island when mutated "No cpg or N"
    #x1 <- n$MeanFreq[YCpG] #create a Value that looks at mean frequency by yes cpg
    #x2 <- n$MeanFreq[NCpG] #create a Value that looks at mean frequency by no cpg
    plot.default( n$MeanFreq + smallestvalue/5, 
        #y values not listed this is intended as x value to be plotted represents two sets of values on a log.
                 xlab = "Position Location", ylab = "Frequency", main = "Location vs frequency CpG non-CpG Graph", # labled titles for plot
                 col = (as.integer(n$makesCpG)+1), #colors the values by their integer, using default association. 
                 log = "y", #the plot is put into a log scale 
                yaxt="n"
    )#close plot.default
    axis(side = 2, at = smallestvalue/5+c(0, round(smallestvalue,4),max(n$MeanFreq)), labels = c(0, round(smallestvalue,4),max(n$MeanFreq)))
    legend('topright', #add legend in top right corner
           legend = c(unique(n$makesCpG)), #names of each category based on the factors in n$makesCpG
           col=unique(as.integer(n$makesCpG)), #colors the factors as is found by integers in n$MakeCpG
           pch=21 #symbols matching dataframe's factors
    ) #close legend
    
    return()# prints out plot with legend 
}#close function

LvsF_CpG_Printer(n) #run function

