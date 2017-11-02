#script for Team NTROPY's function graphing Frequency of synonymous/nonsynonmous/other mutation vs location

#set directory and open data file
setwd("C:/Users/alber/Downloads")

#read csv file into variable "data"
data<-read.csv("OverviewSelCoeff_BachelerFilter.csv")

#plot data mut freq vs location. 
#   Assumptions: plot based on column names, "num", "MeanFreq", "TypeOfSite".
#   Your data must be in these columns! 
#   x axis is location of nucleotide
#   y axis is log of mean frequency
# will have error message to ignore 0's in data because y axis will be translated as logrithmic scale
plotmeanfreqloc<-function(DF){
    if (which(names(DF)=="MeanFreq")==0){
        print("No MeanFreq column! Meanfreq has to be a number")
        return(0)}
    if (which(names(DF)=="Num")==0){
        print("No Num column! Num has to be an integer")}
    if (which(names(DF)=="TypeOfSite")==0){
        print("No TypeOfSite column! Has to be a factor ")}
    
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
    grid(nx = NA, ny = NULL, col = "black", lty = "dotted",
         lwd = par("lwd"), equilogs = TRUE),
    #supress y axis drawing by plot fxn, put # in front to not supress
    yaxt="n",
    # or do log of y axis, delete # symbol
    log="y"
)

#axis function to write y axis with only scale by 10s. 
#USE YOUR OWN APPROPRIATE NUMBERS
axis(
    #2 is to specify left axis (aka y axis)
    2,
    at=c(0.0001,0.001,0.01,0.1)
    # dont need to define labels if same as numbers on axis
    #labels=aty
    #tck marks
    #tck=-0.01
)

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
       text.width=10,
       #justify text legend, xjust=0 is left justified, xjust=1 means right justified
       xjust=1
)
}
