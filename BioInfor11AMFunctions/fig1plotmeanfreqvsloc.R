#script for Team NTROPY's function graphing Frequency of synonymous/nonsynonmous/other mutation vs location

#plot data mut freq vs location and legend
#   Assumptions: plot based on column names, "num", "MeanFreq", "TypeOfSite".
#   Your data must be in these columns! 
#   x axis is location of nucleotide
#   y axis is log of mean frequency
# will have error message to ignore 0's in data because y axis will be translated as logrithmic scale
plotmeanfreqloc<-function(DF, Title){
    
    plot(
        #x vector
        DF$num,
        #y vector
        DF$MeanFreq,
        #make black empty circles as symbol
        pch=21,
        #make outline of symbol black
        col= "black",
        #fill inside of point with color by factor category "TypeOfSite" bg=
        bg=DF$TypeOfSite,
        #Title label
        main = Title,
        #x axis label
        xlab ="Nucleotide Location on Sequence", 
        #y axis label
        ylab ="Mean Frequency of Mutation",
        #cex changes point size
        cex=2,
        #grid superimposes grid onto plot 
        #nx and ny describes x and y axis 
        #NA will automatically set x-axis to default plot ticks 
        grid(nx = NA, ny = NULL, col = "black", lty = "dotted",
             lwd = par("lwd"), equilogs = TRUE),
        #supress y axis drawing by plot fxn, put # in front to not supress
        #yaxt="n"
        # or do log of y axis, delete # symbol
        log="y"
    )
    
    #axis function to write y axis with only scale by 10s. 
    #USE YOUR OWN APPROPRIATE NUMBERS
    #axis(2, at=c(0.0001,0.001,0.01,0.1) , labels=aty, tck=-0.01)
    
    #add legend in top right corner
    legend("topright", 
           #inset legend off from border
           inset= 0.01,
           #names of each category based on factors, alphabetical order of category 1-5 of TypeOfSites
           legend = levels(DF$TypeOfSite), 
           #symbols matching dataframe's factors 1:5 of DF$TypeOfSite
           pch=21,
           #colors matching dataframe's factors 1:5 of DF$TypeOfSite
           col="black",
           #fill colors of circle matching points of plot
           pt.bg=c(1:5),
           #specify scale of whole box of legend to not block data
           cex=.75,
           #specify point size in legend
           pt.cex=3,
           #remove legend border if "n" is specificed. "o" displays border
           bty="o",
           #specific box border thickness/width
           box.lwd=2,
           #specify box border type
           #box.lty=,
           #text.width change
           #text.width=10,
           #justify text legend, xjust=0 is left justified, xjust=1 means right justified
           xjust=1
    )
}
#run fuction
plotmeanfreqloc(BoNS1df, "Bocavirus Mean Frequeny vs Location")