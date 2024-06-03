install.packages("readr")
library(readr)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("pwalign")
BiocManager::install("Biostrings")
BiocManager::install("seqinr")
install.packages("stringr")



library(stringr)
library(pwalign)
library(Biostrings)
library("seqinr")
library(dplyr)



##### read files 1,3
file1<-readr::read_tsv("D:\\R\\SNP task\\file_1.tsv")
file3<-read.csv("D:\\R\\SNP task\\file_3.csv")


##### clean file 1 from na 
file1<-file1[!is.na(file1[["Query_consensus"]]), ]



##### Read the FASTA file2 convert to df
fasta_file <- readLines("D:\\R\\SNP task\\file_2.fasta")

ids <- c()
sequences <- c()

current_id <- ""
current_sequence <- ""

for (line in fasta_file) {
  if (substr(line, 1, 1) == ">") {
    if (current_id != "") {
      ids <- c(ids, current_id)
      sequences <- c(sequences, current_sequence)
    }
    current_id <- gsub(">", "", line)
    current_id <- gsub("^.*\\|[^|]+\\|[^|]+\\|([^|]+)\\|.*$", "\\1", current_id)
    current_sequence <- ""
  } else {
    current_sequence <- paste(current_sequence, line, sep = "")
    
  }
}

if (current_id != "") {
  ids <- c(ids, current_id)
  sequences <- c(sequences, current_sequence)
}

file2 <- data.frame(ID = ids, Sequence = sequences)


##### positions
# 
# indices_list <- list()
# 
# for (i in seq_len(nrow(file2))) {
#   sequence <- file2$Sequence[i]
#   
#   indices <- c()
#   
#   for (j in seq_len(nchar(sequence))) {
#     if (substr(sequence, j, j) == "_") {
#       indices <- c(indices, j)
#     }
#   }
#   
#   indices_list[[i]] <- indices
#   file2$Sequence[i] <- gsub("_", "", sequence)
# }
# 
# file2$Underscore_Positions <- indices_list
# 






##### read file 4 and convert to df
meme_file <- "D:/R/SNP task/file_4.meme"

meme_lines <- readLines(meme_file)
parse_meme <- function(meme_lines) {
  motifs <- list()
  motif_name <- ""
  motif_id <- ""
  matrix_started <- FALSE
  matrix <- NULL
  
  for (line in meme_lines) {
    if (startsWith(line, "MOTIF")) {
      parts <- strsplit(line, " ")[[1]]
      motif_name <- parts[2]
      motif_id <- parts[3]
    } else if (startsWith(line, "letter-probability matrix")) {
      matrix_started <- TRUE
      matrix <- matrix(ncol = 4, nrow = 0)  
    } else if (matrix_started) {
      if (nchar(trimws(line)) == 0) {
        matrix_started <- FALSE
        motifs[[motif_id]] <- list(name = motif_name, matrix = matrix)
        motif_name <- ""
        motif_id <- ""
        matrix <- NULL
      } else {
        row <- as.numeric(unlist(strsplit(trimws(line), "\\s+")))
        if (length(row) == 4) { 
          matrix <- rbind(matrix, row)
        }
      }
    }
  }
  
  if (!is.null(matrix)) {
    motifs[[motif_id]] <- list(name = motif_name, matrix = matrix)
  }
  
  return(motifs)
}

motifs <- parse_meme(meme_lines)



motif_df <- do.call(rbind, lapply(names(motifs), function(motif_id) {
  motif <- motifs[[motif_id]]
  matrix <- motif$matrix
  if (is.null(matrix)) {
    return(NULL)
  }
  data.frame(
    motif_id = motif_id,
    motif_name = motif$name,
    position = 1:nrow(matrix),
    A = matrix[,1],
    C = matrix[,2],
    G = matrix[,3],
    T = matrix[,4],
    stringsAsFactors = FALSE
  )
}))






##### Reverse negative orientation

file1$Target_consensus <- sapply(1:nrow(file1), function(i) {
  if (file1$Orientation[i] == "-") {
    as.character(reverseComplement(DNAString(file1$Target_consensus[i])))
  } else {
    file1$Target_consensus[i]
  }
})







##### overlapping


alignment1 <- data.frame(
  ref_ID = integer(),
  Query_ID = integer(),
  Pattern_Start = integer(),
  Pattern_End = integer(),
  Subject_Start = integer(),
  Subject_End = integer(),
  Alignment_Length = integer(),
  seq = character(),
  stringsAsFactors = FALSE
)

  for (j in seq_len(nrow(file2))) {
    ref <- DNAString(file2$Sequence[j])
    query<-DNAString(file1$Query_ID[1])
    
    alignment <- pairwiseAlignment(ref, query, type = "local")
    
    pattern<-pattern(alignment)
     pattern_start <- start(pattern)
     alignment_length <- width(pattern)
     pattern_end <- pattern_start + alignment_length - 1  
     subject_start <- start(subject(alignment))
     subject_end <- subject_start + alignment_length - 1  
     seq <- as.character(subject(alignment))
     
     
    if (101 >= pattern_start & 101<= pattern_end) {
    
      alignment1 <- rbind(alignment1, data.frame(
         ref_ID = file2$ID[j],
         Query_ID = file1$Query_ID[1],
         Pattern_Start = pattern_start,
         Pattern_End = pattern_end,
         Subject_Start = subject_start,
         Subject_End = subject_end,
         seq = seq,
         Alignment_Length = alignment_length,
         stringsAsFactors = FALSE
       ))
    }
     }




alignment2 <- data.frame(
  target_ID = integer(),
  Query_ID = integer(),
  Pattern_Start = integer(),
  Pattern_End = integer(),
  Subject_Start = integer(),
  Subject_End = integer(),
  Alignment_Length = integer(),
  seq = character(),
  score = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(file1))) {
  Query <- DNAString(file1$Query_ID[1])
  target <- DNAString(file1$Target_consensus[i])
  
  
    
    alignment <- pairwiseAlignment(target, Query, type = "local",gapOpening=10, gapExtension=0.5)
    
     pattern_start <- start(pattern(alignment))
     alignment_length <- width(pattern(alignment))
     pattern_end <- pattern_start + alignment_length - 1  
     subject_start <- start(subject(alignment))
     subject_end <- subject_start + alignment_length - 1  
     seq <- as.character(subject(alignment))
     score<-score(alignment)
    
    if (subject_start <=5 & subject_end >= 5) {
  
     alignment2 <- rbind(alignment2, data.frame(
       target_ID = file1$Target_ID[i],
       Query_ID = file1$Query_ID[i],
       Pattern_Start = pattern_start,
       Pattern_End = pattern_end,
       Subject_Start = subject_start,
       Subject_End = subject_end,
       seq = seq,
       score=score,
       Alignment_Length = alignment_length,
       stringsAsFactors = FALSE
     ))
  }

}



#####  regulatory and non-regulatory

regulatory_df <- data.frame(
  MOTIF_ID = character(),
  regulatory=logical(),
  allele = character(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(alignment2)) {
  positions <- alignment2$Pattern_Start[i] + (5 - alignment2$Subject_Start[i])
  
  
  for (j in 1:nrow(motif_df)) {  
    if (alignment2$target_ID[i] == motif_df$motif_name[j] & motif_df$position[j] == positions) {
      if(motif_df$A[j] == 1){
        regulatory_df <- rbind(regulatory_df, data.frame(
          MOTIF_ID = alignment2$target_ID[i],
          regulatory=TRUE,
          allele = "A",
          stringsAsFactors = FALSE
        ))
      }
       else if(motif_df$C[j]>=1){
         regulatory_df <- rbind(regulatory_df, data.frame(
           MOTIF_ID = alignment2$target_ID[i],
           regulatory=TRUE,
           allele = "C",
           stringsAsFactors = FALSE
         ))
      }
      else if(motif_df$G[j]>=1){
        regulatory_df <- rbind(regulatory_df, data.frame(
          MOTIF_ID = alignment2$target_ID[i],
          regulatory=TRUE,
          allele = "G",
          stringsAsFactors = FALSE
        ))
      }
      else if(motif_df$T[j]>=1){
        regulatory_df <- rbind(regulatory_df, data.frame(
          MOTIF_ID = alignment2$target_ID[i],
          regulatory=TRUE,
          
          allele = "T",
          stringsAsFactors = FALSE
        ))
      }
        else { 
          regulatory_df <- rbind(regulatory_df, data.frame(
            MOTIF_ID = alignment2$target_ID[i],
            regulatory = FALSE,
            allele = " ",
            stringsAsFactors = FALSE
          ))
        }
      }
    }
}



write.csv(regulatory_df, file = "D:\\R\\SNP task\\regulatory_df.csv", row.names = FALSE)
write.csv(alignment1, file = "D:\\R\\SNP task\\alignment1.csv", row.names = FALSE)
write.csv(alignment2, file = "D:\\R\\SNP task\\alignment2.csv", row.names = FALSE)


