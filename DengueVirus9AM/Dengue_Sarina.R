setwd("C:/Users/saqin/Desktop/Class_11AM")
HIVDATA<-read.csv("HIVDATA.csv")
View(HIVDATA)

ColumnG<- HIVDATA[1:985, 7]
head(ColumnG)

}
nucC<-which(ColumnG=="c")
nucA<-which(ColumnG=="a")

column<-data.frame(ColumnG)

paste(ColumnG,collapse = "")->HUH
HUH
gregexpr(pattern = "ca",HUH)->PEN
PEN
PEN->HIVDATA$John
HIVDATA$John<- 0
View(HIVDATA)
PEN
replace(HIVDATA$John(PEN))
View(HIVDATA)


