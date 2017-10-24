#sets the working environment
setwd("C:/Users/saqin/Desktop/Class_11AM")
#reads in csv file and sets to variable HIVDATA
HIVDATA<-read.csv("HIVDATA.csv")
#View the HIV csv file
View(HIVDATA)
#select for column G in dataset and set it to variable ColumnG
ColumnG<- HIVDATA[1:985, 7]


nucC<-which(ColumnG=="c")
nucA<-which(ColumnG=="a")

column<-data.frame(ColumnG)
#collapses data into a string and sets to variable HUH
paste(ColumnG,collapse = "")->HUH
#checks data in HUH
HUH
#looks for pattern ca in the string data HUH, sets information to variable PEN
gregexpr(pattern = "ca",HUH)->PEN
#checks data in PEN
PEN
#makes a new column in HVDATA called John 
PEN->HIVDATA$John
#sets the Column data in John to 0
HIVDATA$John<- 0
View(HIVDATA)


