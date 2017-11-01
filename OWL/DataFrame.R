#Creating a data frame
df<-data.frame(matrix(nrow=2,ncol = 2))
names(df)<-c("hi","bye")
print(df)

newDF<-df
print(newDF)
newDF$hi<-c(1,2)
newDF$bye<-c("Hey")
print(newDF)
