meanFreq <- function(fasta_file){
  library(seqinr)
  dengue_basic <- read.fasta(fasta_file)
  number_of_seqs <- length(dengue_basic)
  dengue_align <- read.alignment(fasta_file, format = "fasta", forceToLower = T)
  dengue_consensus <- seqinr :: consensus(dengue_align, method = "majority")
  dengue_consensus_matrix <- seqinr :: consensus(dengue_align, method = "profile")
  consensus_length <- length(dengue_consensus)
  number_column <- seq(1, consensus_length)
  Dengue_DF <- data.frame("num" = number_column, "MeanFreq" = 0, "WTnt" = dengue_consensus)
  for(x in 1:consensus_length){
    current_base <- dengue_consensus[x]
    current_matrix_base_count <- dengue_consensus_matrix[,x]
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
    Dengue_DF[x, 2] <- ts_count/number_of_seqs
  }
  return(Dengue_DF)
}

#example
#meanFreq("DengueViruses/DengueVirus1.fasta")
