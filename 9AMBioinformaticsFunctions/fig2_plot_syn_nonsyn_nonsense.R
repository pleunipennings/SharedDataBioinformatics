#input should be dataframe
#the command should be something like the function pulls out information
  #from the dataframe and uses that info to plot
#output should be plot of Meanfreq vs. Position (num), colored by TypeofSite
plot_NonSynNon = function(df) {
    plot(df$num, df$freq+0.1, type = "p", main = "Nucleotide Positon (num) vs. Mean Frequency (freq) of Virus", 
         log = "y", xlab = "num", ylab = "log of freq", col = as.factor(df$TypeOfSite))
}


#explanation of some plot arguments
    #df$freq+0.1 we had to do +0.1 because our graph looked better that way, yours may or may not need more zeros
    #type = "p" makes a scatterplot
    #col = as.factor(df$TypeOfSite) was so that we could color everything based on TypeOfSite, but first we had to 
        #ensure that R read it as factors; otherwise it wouldn't color it b/c "nonsyn" was not a color (tl;dr R is picky)


virus_plot(final_BOCA) #this is what using the function should look like, our df was titled final_BOCA


#making a legend
legend("topleft",
       c("syn", "nonsyn", "nonsense"),
       lty = c(1,1), lwd=4,
       col=c("darkolivegreen3", "red", "purple"),
       bty = "n")

#who is in the team: Shannel, Milo, Ricky, Natalie
#who helped make the function: mostly Shannel, sprinkles of Milo and Ricky here and there
