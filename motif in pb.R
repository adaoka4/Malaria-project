library(rtracklayer)
library(Biostrings)
library(purrr)

# Load gff and extract utr features and gene features
gff_file <- "PlasmoDB-63_PbergheiANKA.gff"
gff <- readGFF(gff_file)

gene_features <- subset(gff, type == "protein_coding_gene")


# Load fasta file and convert to a data.frame
fasta_file <- "PlasmoDB-63_PbergheiANKA_Genome.fasta"
fasta_sequences <- rtracklayer::import(fasta_file)
fasta_sequences2 <- as.data.frame(fasta_sequences)

# Extract first text before "|" symbol from row names
new_row_names <- sub("\\|.*", "", rownames(fasta_sequences2))
new_row_names <- sub("\\s*\\|.*", "", rownames(fasta_sequences2))

# Assign new row names to the dataframe and convert rownames to a column
rownames(fasta_sequences2) <- new_row_names
library(tibble)
fasta_sequences2 <- rownames_to_column(fasta_sequences2, var="seqid")


# Convert gene_features to a data.frame and select relevant columns
gene_features_df <- as.data.frame(gene_features)
gene_select <- c(1:5,7,9:12)
gene_features_df <- gene_features_df[, gene_select]


# Function to flip and substitute characters in a sequence
flip_minus_sequence <- function(Sequence, strand) {
  if (strand == "-") {
    reversed_sequence <- paste(rev(strsplit(Sequence, "")[[1]]), collapse = "")
    flipped_sequence <- strsplit(reversed_sequence, "")[[1]]

    for (i in seq_along(flipped_sequence)) {
      letter <- flipped_sequence[i]

      if (letter == "A") {
        flipped_sequence[i] <- "T"
      } else if (letter == "T") {
        flipped_sequence[i] <- "A"
      } else if (letter == "C") {
        flipped_sequence[i] <- "G"
      } else if (letter == "G") {
        flipped_sequence[i] <- "C"
      }
    }

    flipped_sequence <- paste(flipped_sequence, collapse = "")
  } else {
    flipped_sequence <- Sequence
  }

  return(flipped_sequence)
}


# Function to extract gene sequences based on start and end positions
extract_gene_sequences <- function(fasta_sequences2, gene_features_df) {
  sequences <- vector("list", nrow(gene_features_df))
  
  for (i in 1:nrow(gene_features_df)) {
    row <- gene_features_df[i, ]
    seq_id <- row$seqid
    start_pos <- row$start
    end_pos <- row$end
    
    dna_seq <- fasta_sequences2[fasta_sequences2$seqid == seq_id, "x"]
    sequence <- substring(dna_seq, start_pos, end_pos)
    
    row$Sequence <- sequence
    sequences[[i]] <- row
  }
  
  sequences_df <- do.call(rbind, sequences)
  return(sequences_df)
}

# Call the function to extract sequences
gene_sequences <- extract_gene_sequences(fasta_sequences2, gene_features_df)

# Apply the flip_sequence function to the Sequence and Strand columns in gene_sequences
gene_sequences$FlippedSequence <- mapply(flip_minus_sequence, gene_sequences$Sequence, gene_sequences$strand)
gene_sequences$Parent <- purrr::map_chr(gene_sequences$Parent, toString)

# Save dataframe to a CSV file
write.csv(gene_sequences, file = "~/Pb_genes.csv", row.names = FALSE)


# Convert data frames to fasta
library(seqinr)

# genes
# Define the output fasta file path
gene_output_file <- "Pb_genes.fasta"

# Open the output file
gene_file_out <- file(gene_output_file, "w")

# Write each sequence to the output file
for (i in 1:nrow(gene_sequences)) {
  # Write the gene name
  writeLines(paste(">", gene_sequences$ID[i], sep="",
                   " | ", gene_sequences$Name[i], 
                   " | ", gene_sequences$description[i],
                   " | start = ", gene_sequences$start[i],
                   " | end = ", gene_sequences$end[i],
                   " | strand = ", gene_sequences$strand[i]), gene_file_out)
  
  # Write the sequence
  writeLines(gene_sequences$FlippedSequence[i], gene_file_out)
}

# Close the output file
close(gene_file_out)




