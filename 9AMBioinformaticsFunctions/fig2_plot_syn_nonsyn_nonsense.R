#input should be dataframe
#the command should be something like the function pulls out information
#from the dataframe and uses that info to plot
#output should be plot of Meanfreq vs. Position (num), colored by TypeofSite
virus_plot = function(df) {
    plot(df$num, df$freq+0.1, type = "p", main = "Nucleotide Positon (num) vs. Mean Frequency (freq) of Virus", 
         log = "y", xlab = "num", ylab = "log of freq", col = as.factor(df$TypeOfSite))
}

virus_plot(final_BOCA) #THIS FINALLY WORKED???

#making a legend
legend("topleft",
       c("syn", "nonsyn", "nonsense"),
       lty = c(1,1), lwd=4,
       col=c("darkolivegreen3", "red", "purple"),
       bty = "n")