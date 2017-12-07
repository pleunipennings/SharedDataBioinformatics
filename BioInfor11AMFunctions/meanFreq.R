meanFreq <- function(fasta_file){
  library(seqinr)
  virus_basic <- read.fasta(fasta_file)
  number_of_seqs <- length(virus_basic)
  virus_align <- read.alignment(fasta_file, format = "fasta", forceToLower = T)
  virus_consensus <- seqinr :: consensus(virus_align, method = "majority")
  virus_consensus_matrix <- seqinr :: consensus(virus_align, method = "profile")
  consensus_length <- length(virus_consensus)
  number_column <- seq(1, consensus_length)
  virus_DF <- data.frame("num" = number_column, "MeanFreq" = 0, "wtnt" = virus_consensus)
  for(x in 1:consensus_length){
    current_base <- virus_consensus[x]
    current_matrix_base_count <- virus_consensus_matrix[,x]
    ts_count <- 0
        if(current_base == "a"){
      ts_count <- current_matrix_base_count[["g"]]
    }
    if(current_base == "g"){
      ts_count <- current_matrix_base_count[["a"]]
    }
    if(current_base == "c"){
      ts_count <- current_matrix_base_count[["t"]]
    }
    if(current_base == "t"){
      ts_count <- current_matrix_base_count[["c"]]
    }
    virus_DF[x, 2] <- ts_count/number_of_seqs
  }
  virus_DF$wtnt<-as.character(virus_DF$wtnt)
  return(virus_DF)
}

#made by 11 AM dengue team

#example
#meanFreq("virusViruses/virusVirus1.fasta")
