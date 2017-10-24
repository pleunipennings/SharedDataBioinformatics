#script for Team NTROPY's function graphing Frequency of synonymous/nonsynonmous/other mutation vs location

#set directory and open data file
setwd("C:/Users/alber/Downloads")

#read csv file into variable "data"
data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")


#set overall scale of plots/legends
op <- par(cex = 1)

#plot data mut freq vs location. 
#   Assumptions: plot based on column names, "num", "MeanFreq", "TypeOfSite".
#   Your data must be in these columns! 
plot(
    #x vector
    data$num,
    #y vector
     log(data$MeanFreq),
      #make black empty circles as symbol
     pch=21,
    #make outline of symbol black
     col= "black",
    #fill inside of point with color by factor category "TypeOfSite" bg=
    bg=data$TypeOfSite,
     #Title label
    main = "HIV Practice Data",
    #x axis label
    xlab ="Location on Sequence", 
    #y axis label
    ylab ="Mean Frequency of Mutation",
    #cex changes point size
    cex=2,
    #grid superimposes grid onto plot 
    #nx and ny describes x and y axis 
    #NA will automatically set x-axis to default plot ticks 
    grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
)

#set overall scale of plots/legends
op <- par(cex = 1)

#add legend in top right corner
legend("topright", 
       #inset legend off from border
       inset= 0.01,
       #names of each category based on factors, alphabetical order of category 1-5 of TypeOfSites
       legend = levels(data$TypeOfSite), 
       #symbols matching dataframe's factors 1:5 of data$TypeOfSite
       pch=21,
       #colors matching dataframe's factors 1:5 of data$TypeOfSite
       col="black",
       #fill colors of circle matching points of plot
       pt.bg=c(1:5),
       #specify scale of whole box of legend to not block data
       cex=1,
       #specify point size in legend
       pt.cex=2,
       #remove legend border if "n" is specificed. "o" displays border
       bty="o",
       #specific box border thickness/width
       box.lwd=2,
       #specify box border type
       #box.lty=,
       #text.width change
       text.width=30,
       #justify text legend, xjust=0 is left justified, xjust=1 means right justified
       xjust=0.5

)
