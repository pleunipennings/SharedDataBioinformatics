read.csv("HIVDATA.CSV") ->HIV


# Collapses data into a string and sets to variable STRING
paste(HIV$WTnt, collapse = '') -> STRING
STRING

#Looks for pattern tg in the data STRING, sets location sites to variable TG
gregexpr(pattern ='tg',STRING ) -> TG
TG

# Insert LIST TG into a data frame BELL
BELL <- data.frame(matrix(unlist(TG)))
BELL

# create new column name John with values zero
HIV$JOHN<- 0

# Inserting value 1 in every TG site using BELL
HIV[BELL[,1],"JOHN"] <- 1


# For CA sites: 

#Looks for pattern ca in the data STRING, sets location sites to variable CA
gregexpr(pattern ='ca',STRING ) -> CA
CA

# Insert LIST CA into a data frame CASITES
CASITES <- data.frame(matrix(unlist(CA)))
CASITES

# Inserting value 1 in every CA site using CA
HIV[CASITES[,1]+1,"JOHN"] <- 1



