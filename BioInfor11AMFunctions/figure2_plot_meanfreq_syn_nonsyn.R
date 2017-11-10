#11AM function to plot mean freq vs location with synoymous and non synonmous Type of Sites
#   Team: Liana, Alejandro, Albert

#HOW TO USE: 
#  plotsyn(DF, "yourtitlehere", logy)
#   1st argument, DF is your dataframe 
#   2nd argument, Insert the title you want in "Title." Make sure to use quotes!
#       e.g. plotmeanfreqloc(DF, "your title here")
#   3rd argument, type logy no quotes if you want y axis to be logrithmic
#       e.g. plotmeanfreqloc(dataframe, "your graph title", logy)
#There will be an error message about infinite axis limits but it will still work.
#ASSUMPTIONS: requires dataframe column names, "num", "MeanFreq", "TypeOfSite".
#       num and MeanFreq column have to be numbers
#        TypeOfSite has to be factors
#NOTE: Code will omit 0's in MeanFreq. Presumably, the 0's do not have 
#   biological significance for this figure

plotsyn<-function(DF, Title, logy){
    #code to ignore any 0's in MeanFreq column
    DF->df2
    #df1<-data.frame("num"=df$num, "MeanFreq"=df$MeanFreq, "TypeOfSite"=df$TypeOfSite)
    #NA->df1[df1==0]
   # na.omit(df1)->df2
#set xy window to fit data from 0 to its max
        plot.window(c(0,nrow(DF)), c(0,max(DF$MeanFreq)))
    if(missing("logy")){
        plot(
             #x vector
            df2$num,
            #y vector
            df2$MeanFreq,
            #make symbols to each factor of TypeOfSite
            pch=c(df2$TypeOfSite),
            #make symbol colors to each factor of Type Of Site
            col=df2$TypeOfSite,
            #Title label
            main = Title,
            #x axis label
            xlab ="Nucleotide Location on Sequence", 
            #y axis label
            ylab ="Mean Frequency of Mutation",
            #cex changes overall magnification
            cex=2,
            #grid superimposes gridlines onto plot 
            #nx and ny describes x and y axis 
            #NA will automatically set x-axis to default plot ticks 
            grid(nx = NA, ny = NULL, col = "black", lty = "dotted",
                 lwd = par("lwd"), equilogs = TRUE)
         )
    }
    else {
        plot(
            #x vector
            df2$num,
            #y vector
            df2$MeanFreq,
            #make symbols to each factor of TypeOfSite
            pch=c(df2$TypeOfSite),
            #make symbol colors to each factor of Type Of Site
            col=df2$TypeOfSite,
            #Title label
            main = Title,
            #x axis label
            xlab ="Nucleotide Location on Sequence", 
            #y axis label
            ylab ="Mean Frequency of Mutation",
            #cex changes overall magnification
            cex=2,
            #grid superimposes gridlines onto plot 
            #nx and ny describes x and y axis 
            #NA will automatically set x-axis to default plot ticks 
            grid(nx = NA, ny = NULL, col = "black", lty = "dotted",
                 lwd = par("lwd"), equilogs = TRUE),
            log="y"
        )
    }
    
    #add legend in top right corner
    legend("topright", 
           #inset legend off from border
           inset= 0.01,
           #names of each category based on factors, alphabetical order of category 1-5 of TypeOfSites
           legend = levels(df2$TypeOfSite), 
           #symbols matching dataframe's factors 1:5 of DF$TypeOfSite
           pch=c(1:5),
           #colors matching dataframe's factors 1:5 of DF$TypeOfSite
           col=c(1:5),
           #specify scale of whole box of legend to not block data
           cex=.75,
           #specify point size in legend
           pt.cex=3,
           #remove legend border if "n" is specificed. "o" displays border
           bty="o",
           #specific box border thickness/width
           box.lwd=2,
           #text.width change
           text.width=120,
           #justify text legend, xjust=0 is left justified, xjust=1 means right justified
           xjust=0
    )
}

#   Team: Liana, Alejandro, Albert
