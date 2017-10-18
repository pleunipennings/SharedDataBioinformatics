# freq calc v2 GL

library(seqinr)

bk <- read.fasta("bk.txt")

# reference sequence
ref <- bk[[1]]

# dataframe columns
num <- c(1:1089)
wtnt <- c()
freq <- c()

# for freq calculation later
absfreq <- c(rep(0,1089))
totalcount <- c(rep(0,1089))

# average WT calculation
# counts number of each nucleotide in each position
acount <- c(rep(0,1089))
gcount <- c(rep(0,1089))
ccount <- c(rep(0,1089))
tcount <- c(rep(0,1089))
nuc <- c()

# same as line 20 comment
for (i in 1:length(bk)) {
  sequence <- bk[[i]]
  for (j in 1:length(sequence)) {
    if (sequence[j] == 'a') {acount[j] = acount[j] + 1}
    if (sequence[j] == 'g') {gcount[j] = gcount[j] + 1}
    if (sequence[j] == 'c') {ccount[j] = ccount[j] + 1}
    if (sequence[j] == 't') {tcount[j] = tcount[j] + 1}
  }
}

# assigns wtnt based on most frequent nucleotide across all sequences per position
for (j in 1:length(sequence)) {
  nuc[j] <- max(c(acount[j],gcount[j],ccount[j],tcount[j]))
  if (max(nuc[j]) == acount[j]) {wtnt[j] <- 'a'}
  if (max(nuc[j]) == gcount[j]) {wtnt[j] <- 'g'}
  if (max(nuc[j]) == ccount[j]) {wtnt[j] <- 'c'}
  if (max(nuc[j]) == tcount[j]) {wtnt[j] <- 't'}
  nuc <- c()
}

# gives absolute totals to be used for frequency calculation
for (i in 1:length(bk)) {
  sequence <- bk[[i]]
  for (j in 1:length(sequence)){
    if (wtnt[j] == 'a') {
      if (sequence[j] == 'g') {
        absfreq[j] <- absfreq[j] + 1
        totalcount[j] <- totalcount[j] + 1
      }
      else if (sequence[j] == 'a') {
        totalcount[j] <- totalcount[j] + 1
      }
    }
    if (wtnt[j] == 'g') {
      if (sequence[j] == 'a') {
        absfreq[j] <- absfreq[j] + 1
        totalcount[j] <- totalcount[j] + 1
      }
      else if (sequence[j] == 'g') {
        totalcount[j] <- totalcount[j] + 1
      }
    }
    if (wtnt[j] == 'c') {
      if (sequence[j] == 't') {
        absfreq[j] <- absfreq[j] + 1
        totalcount[j] <- totalcount[j] + 1
      }
      else if (sequence[j] == 'c') {
        totalcount[j] <- totalcount[j] + 1
      }
    }
    if (wtnt[j] == 't') {
      if (sequence[j] == 'c') {
        absfreq[j] <- absfreq[j] + 1
        totalcount[j] <- totalcount[j] + 1
      }
      else if (sequence[j] == 't') {
        totalcount[j] <- totalcount[j] + 1
      }
    }
  }
}

# calculates frequency as percentage
for (i in 1:length(absfreq)) {
  freq[i] <- absfreq[i] / totalcount[i]
}

# creates dataframe containing all data
bk_data <- data.frame(num,wtnt,freq)

# outputs data to file
write.csv(bk_data,"bk_data.csv")
plot(bk_data$freq + 0.01, log="y", col=c(1,2,3))
#1 black fist aa in codon, 2 red 2nd aa in codon, 3 green 3rd aa in codon
