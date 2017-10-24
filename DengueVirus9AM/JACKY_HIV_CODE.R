read.csv("HIVDATA.CSV") ->HIV

HIV

HIV$WTnt

paste(HIV$WTnt, collapse = '') -> STRING

STRING


gregexpr(pattern ='tg',STRING ) -> TG
TG

BELL <- data.frame(matrix(unlist(TG)))
BELL

HIV$NUM1 = HIV$num
HIV$WTnt1 = HIV$WTnt

HIV$JOHN<- 0

HIV[BELL[,1],"JOHN"] <- 1


# For CA sites: 

gregexpr(pattern ='ca',STRING ) -> CA
CA

CASITES <- data.frame(matrix(unlist(CA)))
CASITES


HIV[CASITES[,1]+1,"JOHN"] <- 1



