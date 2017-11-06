
## Group # 1 ( 9am) - CPG Sites Function
## Contributior - Jacky, Sarina, Preeti


######## Investigating CPG sites: 

## Note : this function assumes that wtnt column is present in the 
##        dataframe and your dataframe is set as df

# Runs CPG_site function which takes in the arguement/input dataframe
CPG_site<-function(df){
# Collapses wtnt vector sequence into a character string and sets to variable STRING
paste(df$wtnt, collapse = '') -> STRING
STRING

#Looks for pattern tg in the data STRING, return location of TG sites within the STRING and store in variable TG
gregexpr(pattern ='tg',STRING ) -> TG
TG

# Since TG is a list, insert 'list' TG into a data frame BELL
BELL <- data.frame(matrix(unlist(TG)))
BELL

# create new column name CPG with values zero
df$CPG<- 0

# Inserting value 1 in the column "CPG" and using the values found dataframe BELL as the row #
df[BELL[,1],"CPG"] <- 1


# For CA sites: 

#Looks for pattern ca in the data STRING, return location of CA sites within the STRING and store in variable CA
gregexpr(pattern ='ca',STRING ) -> CA
CA

# Since TG is a list, insert 'list' TG into a data frame BELL
CASITES <- data.frame(matrix(unlist(CA)))
CASITES

# Inserting value 1 in the column "CPG" and using the values found dataframe BELL as the row #
df[CASITES[,1]+1,"CPG"] <- 1
}

# Calls the CPG_site function with the dataframe as the argument
CPG_site(df)




