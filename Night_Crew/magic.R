# freq calc v2 GL

library(seqinr)
library(stringi)

bk <- read.fasta("bk.txt")
bk<-read.fasta("InfluenzaAvirus_HA_H1N1.txt")
nightcrewBK= function(data) {
# dataframe columns
num <- c(1:1089)
WTnt <- c()
MeanFreq <- c()

# for MeanFreq calculation later
absfreq <- c(rep(0, 1089))
totalcount <- c(rep(0, 1089))

# average WT calculation
# counts number of each nucleotide in each position
acount <- c(rep(0, 1089))
gcount <- c(rep(0, 1089))
ccount <- c(rep(0, 1089))
tcount <- c(rep(0, 1089))
nuc <- c()

# same as line 20 comment
for (i in 1:length(bk)) {
    sequence <- bk[[i]]
    for (j in 1:length(sequence)) {
        if (sequence[j] == 'a') {
            acount[j] = acount[j] + 1
        }
        if (sequence[j] == 'g') {
            gcount[j] = gcount[j] + 1
        }
        if (sequence[j] == 'c') {
            ccount[j] = ccount[j] + 1
        }
        if (sequence[j] == 't') {
            tcount[j] = tcount[j] + 1
        }
    }
}

# assigns wtnt based on most frequent nucleotide across all sequences per position
for (j in 1:length(sequence)) {
    nuc[j] <- max(c(acount[j], gcount[j], ccount[j], tcount[j]))
    if (max(nuc[j]) == acount[j]) {
        WTnt[j] <- 'a'
    }
    if (max(nuc[j]) == gcount[j]) {
        WTnt[j] <- 'g'
    }
    if (max(nuc[j]) == ccount[j]) {
        WTnt[j] <- 'c'
    }
    if (max(nuc[j]) == tcount[j]) {
        WTnt[j] <- 't'
    }
    nuc <- c()
}

# gives absolute totals to be used for frequency calculation
for (i in 1:length(bk)) {
    sequence <- bk[[i]]
    for (j in 1:length(sequence)) {
        if (WTnt[j] == 'a') {
            if (sequence[j] == 'g') {
                absfreq[j] <- absfreq[j] + 1
                totalcount[j] <- totalcount[j] + 1
            }
            else if (sequence[j] == 'a') {
                totalcount[j] <- totalcount[j] + 1
            }
        }
        if (WTnt[j] == 'g') {
            if (sequence[j] == 'a') {
                absfreq[j] <- absfreq[j] + 1
                totalcount[j] <- totalcount[j] + 1
            }
            else if (sequence[j] == 'g') {
                totalcount[j] <- totalcount[j] + 1
            }
        }
        if (WTnt[j] == 'c') {
            if (sequence[j] == 't') {
                absfreq[j] <- absfreq[j] + 1
                totalcount[j] <- totalcount[j] + 1
            }
            else if (sequence[j] == 'c') {
                totalcount[j] <- totalcount[j] + 1
            }
        }
        if (WTnt[j] == 't') {
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
    MeanFreq[i] <- absfreq[i] / totalcount[i]
}
# translation and comparison setup
TypeOfSite <- c()
MUTAA <- c()
WTAAp <- translate(WTnt, NAstring = "X")

#fig out true WTAA
WTAAs <- stri_dup(WTAAp, 3)
WTAA <- unlist(strsplit(WTAAs, ""))

WTAA <- as.character(WTAA)
# creates dataframe containing all data
bk_data <- data.frame(num, WTnt, MeanFreq, WTAA)

}

bk_data<-nightcrewBK(bk)



#syn <- c(rep(0, 1089))
#nonsyn <- c(rep(0, 1089))
#overlap <- c(rep(0, 1089))
#nonsense <- c(rep(0, 1089))
# finds type of site based on average type
# for (i in 1:length(bk)) {
#     seq <- bk[[i]]
#     mutated <- translate(seq, NAstring = "X")
#     for (j in 1:length(seq)) {
#         if (WTnt[j] != seq[j] &&
#             WTAA[((j - 1) %/% 3) + 1] != mutated[((j - 1) %/% 3) + 1]) {
#             if (mutated[((j - 1) %/% 3) + 1] != "*") {
#                 nonsyn[j] <- nonsyn[j] + 1
#             } else {
#                 nonsense[j] <- nonsense[j] + 1
#             }
#         } else if (WTnt[j] != seq[j] &&
#                    WTAA[((j - 1) %/% 3) + 1] == mutated[((j - 1) %/% 3) + 1]) {
#             if (mutated[((j - 1) %/% 3) + 1] != "*") {
#                 syn[j] <- syn[j] + 1
#             } else {
#                 nonsense[j] <- nonsense[j] + 1
#             }
#         } else if (WTnt[j] == seq[j]) {
#             if (mutated[((j - 1) %/% 3) + 1] != "*") {
#                 overlap[j] <- overlap[j] + 1
#             } else {
#                 nonsense[j] <- nonsense[j] + 1
#             }
#         }
#     }
# }
# 
# # determines typeofsite based on most appearances of type
# # problem: skews towards overlap
# type <- 0
# for (j in 1:length(sequence)) {
#     type[j] <- max(c(nonsyn[j], syn[j], overlap[j], nonsense[j]))
#     if (max(type[j]) == nonsyn[j]) {
#         TypeOfSite[j] <- 'nonsyn'
#     }
#     if (max(type[j]) == syn[j]) {
#         TypeOfSite[j] <- 'syn'
#     }
#     if (max(type[j]) == overlap[j]) {
#         TypeOfSite[j] <- 'overlap'
#     }
#     if (max(type[j]) == nonsense[j]) {
#         TypeOfSite[j] <- 'nonsense'
#     }
#     nuc <- c()
# }
# 
# 

# creates dataframe containing all data
#bk_data <- data.frame(num, WTnt, MeanFreq, WTAA)



#finding the MUTAA
# nightcrewMUTAA = function(data) {
#     bk_data$MUTAA <- c(0)
#     #MUTAA
#     x = 0
#     a = 1
#     bk_data$WTnt -> bk_data$A
#     for (i in 1:363) {
#         x <- a
#         #print(x)
#        # print(a)
#         #print(i)
#         if (bk_data$WTnt[x] == "a") {
#             bk_data$A[x] <- "g"
#             a = 3 + a
#         }
#         if (bk_data$WTnt[x] == "g") {
#             bk_data$A[x] <- "a"
#             a = 3 + a
#         }
#         if (bk_data$WTnt[x] == "c") {
#             bk_data$A[x] <- "t"
#             a = 3 + a
#         }
#         if (bk_data$WTnt[x] == "t") {
#             bk_data$A[x] <- "c"
#             a = 3 + a
#         }
#     }
#     
#     x = 0
#     a = 2
#     bk_data$WTnt -> bk_data$B
#     for (i in 1:363) {
#         x <- a
#         #print(x)
#         #print(a)
#         if (bk_data$WTnt[x] == "a") {
#             bk_data$B[x] <- "g"
#             a = 3 + a
#         }
#         if (bk_data$WTnt[x] == "g") {
#             bk_data$B[x] <- "a"
#             a = 3 + a
#         }
#         if (bk_data$WTnt[x] == "c") {
#             bk_data$B[x] <- "t"
#             a = 3 + a
#         }
#         if (bk_data$WTnt[x] == "t") {
#             bk_data$B[x] <- "c"
#             a = 3 + a
#         }
#     }
#     x = 0
#     a = 3
#     bk_data$WTnt -> bk_data$C
#     for (i in 1:363) {
#         x <- a
#         #print(x)
#         #print(a)
#         if (bk_data$WTnt[x] == "a") {
#             bk_data$C[x] <- "g"
#             a = 3 + a
#         }
#         if (bk_data$WTnt[x] == "g") {
#             bk_data$C[x] <- "a"
#             a = 3 + a
#         }
#         if (bk_data$WTnt[x] == "c") {
#             bk_data$C[x] <- "t"
#             a = 3 + a
#         }
#         if (bk_data$WTnt[x] == "t") {
#             bk_data$C[x] <- "c"
#             a = 3 + a
#         }
#     }
#     
#     
#     
#     
#     MUTAA1 <- translate(as.character(bk_data$A), NAstring = "X")
#     MUTAA1 <- stri_dup(MUTAA1, 3)
#     MUTAA1 <- unlist(strsplit(MUTAA1, ""))
#     
#     MUTAA2 <- translate(as.character(bk_data$B), NAstring = "X")
#     MUTAA2 <- stri_dup(MUTAA2, 3)
#     MUTAA2 <- unlist(strsplit(MUTAA2, ""))
#     
#     MUTAA3 <- translate(as.character(bk_data$C), NAstring = "X")
#     MUTAA3 <- stri_dup(MUTAA3, 3)
#     MUTAA3 <- unlist(strsplit(MUTAA3, ""))
#     
#     bk_data$MUTAA1 <- MUTAA1
#     bk_data$MUTAA2 <- MUTAA2
#     bk_data$MUTAA3 <- MUTAA3
#     w = 1
#     e = 2
#     t = 3
#     bk_data$MUTAA = c(0)
#     for (i in 1:363) {
#         v <- w
#         if (bk_data$WTnt[v] == "a" && bk_data$A[v] == "g") {
#             bk_data$MUTAA[v] = as.character(bk_data$MUTAA1[v])
#             w = 3 + w
#         }
#         if (bk_data$WTnt[v] == "g" && bk_data$A[v] == "a") {
#             bk_data$MUTAA[v] = as.character(bk_data$MUTAA1[v])
#             w = 3 + w
#         }
#         if (bk_data$WTnt[v] == "t" && bk_data$A[v] == "c") {
#             bk_data$MUTAA[v] = as.character(bk_data$MUTAA1[v])
#             w = 3 + w
#         }
#         if (bk_data$WTnt[v] == "c" && bk_data$A[v] == "t") {
#             bk_data$MUTAA[v] = as.character(bk_data$MUTAA1[v])
#             w = 3 + w
#         }
#     }
#     for (i in 1:363) {
#         r <- e
#         
#         if (bk_data$WTnt[r] == "a" && bk_data$B[r] == "g") {
#             bk_data$MUTAA[r] = as.character(bk_data$MUTAA2[r])
#             e = 3 + e
#         }
#         if (bk_data$WTnt[r] == "g" && bk_data$B[r] == "a") {
#             bk_data$MUTAA[r] = as.character(bk_data$MUTAA2[r])
#             e = 3 + e
#         }
#         if (bk_data$WTnt[r] == "t" && bk_data$B[r] == "c") {
#             bk_data$MUTAA[r] = as.character(bk_data$MUTAA2[r])
#             e = 3 + e
#         }
#         if (bk_data$WTnt[r] == "c" && bk_data$B[r] == "t") {
#             bk_data$MUTAA[r] = as.character(bk_data$MUTAA2[r])
#             e = 3 + e
#         }
#     }
#     for (i in 1:363) {
#         c <- t
#         if (bk_data$WTnt[c] == "a" && bk_data$C[c] == "g") {
#             bk_data$MUTAA[c] = as.character(bk_data$MUTAA3[c])
#             t = 3 + t
#         }
#         if (bk_data$WTnt[c] == "g" && bk_data$C[c] == "a") {
#             bk_data$MUTAA[c] = as.character(bk_data$MUTAA3[c])
#             t = 3 + t
#         }
#         if (bk_data$WTnt[c] == "t" && bk_data$C[c] == "c") {
#             bk_data$MUTAA[c] = as.character(bk_data$MUTAA3[c])
#             t = 3 + t
#         }
#         if (bk_data$WTnt[c] == "c" && bk_data$C[c] == "t") {
#             bk_data$MUTAA[c] = as.character(bk_data$MUTAA3[c])
#             t = 3 + t
#         }
#     }
#     
#     bk_data <-subset(bk_data, select = -c(MUTAA1, MUTAA2, MUTAA3, A, B, C))
#     return(bk_data)
# }
# bk_data<-nightcrewMUTAA(bk_data)

#simple version of MUTAA
MUTTA= function(bk_data){
x=1
for (x in 1:nrow(bk_data)) {
    
    #print(x)
    # print(a)
    #print(i)
    if (bk_data$WTnt[x] == "a") {
        bk_data$A[x] <- "g"
    }
    if (bk_data$WTnt[x] == "g") {
        bk_data$A[x] <- "a"
    }
    if (bk_data$WTnt[x] == "c") {
        bk_data$A[x] <- "t"
    }
    if (bk_data$WTnt[x] == "t") {
        bk_data$A[x] <- "c"
    }
}

x=1
y=1
count<-1
bk_data$MA<-c(0)
for (x in 1:(nrow(bk_data)/3)) {
    #print(x)
    #print(y)
    for(y in 1:3){
        if(y==1){
            bk_data$MA[count]<-translate(r<-c(bk_data$A[count],as.character(bk_data$WTnt[count+1]),as.character(bk_data$WTnt[count+2])))
        }
        if(y==2){
            bk_data$MA[count]<-translate(a<-c(as.character(bk_data$WTnt[count-1]),as.character(bk_data$A[count]),as.character(bk_data$WTnt[count+1])))
        }
        if(y==3){
            bk_data$MA[count]<-translate(p<-c(as.character(bk_data$WTnt[count-2]),as.character(bk_data$WTnt[count-1]),as.character(bk_data$A[count])))
        }
        count<-count+1
       
    }
    
    }


}

    # outputs data to file
    write.csv(bk_data, "bk_data.csv")
    save(bk_data, file = "bk_data.Rda")
    load("bk_data.Rda")
    
    #check plot
    plot(bk_data$MeanFreq + 0.01,
         log = "y",
         col = c(1, 2, 3))
    #1 black fist aa in codon, 2 red 2nd aa in codon, 3 green 3rd aa in codon
    