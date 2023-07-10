library(rtracklayer)
library(Biostrings)
library(purrr)

# Load gff and extract utr features and gene features
gff_file <- "PlasmoDB-63_Pfalciparum3D7.gff"
gff <- readGFF(gff_file)
utr_features <- subset(gff, type == "five_prime_UTR")
gene_features <- subset(gff, type == "protein_coding_gene")



# Load fasta file and convert to a data.frame
fasta_file <- "PlasmoDB-63_Pfalciparum3D7_Genome.fasta"
fasta_sequences <- rtracklayer::import(fasta_file)
fasta_sequences2 <- as.data.frame(fasta_sequences)

# Extract first text before "|" symbol from row names
new_row_names <- sub("\\|.*", "", rownames(fasta_sequences2))
new_row_names <- sub("\\s*\\|.*", "", rownames(fasta_sequences2))

# Assign new row names to the dataframe and convert rownames to a column
rownames(fasta_sequences2) <- new_row_names
library(tibble)
fasta_sequences2 <- rownames_to_column(fasta_sequences2, var="seqid")

# Convert utr_features to a data.frame and select relevant columns
utr_features_df <- as.data.frame(utr_features)
select <- c(1:5,7,9,13)
utr_features_df <- utr_features_df[, select]


# Convert gene_features to a data.frame and select relevant columns
gene_features_df <- as.data.frame(gene_features)
gene_select <- c(1:5,7,9:12)
gene_features_df <- gene_features_df[, gene_select]

# Function to extract utr sequences based on start and end positions
extract_utr_sequences <- function(fasta_sequences2, utr_features_df) {
  sequences <- vector("list", nrow(utr_features_df))
  
  for (i in 1:nrow(utr_features_df)) {
    row <- utr_features_df[i, ]
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
utr_sequences <- extract_utr_sequences(fasta_sequences2, utr_features_df)

utr_sequences$Parent <- purrr::map_chr(utr_sequences$Parent, toString)

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


# Apply the flip_sequence function to the Sequence and Strand columns in utr_sequences
utr_sequences$FlippedSequence <- mapply(flip_minus_sequence, utr_sequences$Sequence, utr_sequences$strand)
utr_sequences$length <- (utr_sequences$end - utr_sequences$start) + 1

# Save dataframe to a CSV file
#write.csv(utr_sequences, file = "~/Pf_annotated_utr.csv", row.names = FALSE)


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




# Convert data frames to fasta
library(seqinr)

# Sort the dataframe by decreasing lengths
utr_sequences <- utr_sequences[order(-utr_sequences$length), ]

# Define the output fasta file path
utr_output_file <- "Pf_annotated_utr.fasta"

# Open the output file
utr_file_out <- file(utr_output_file, "w")

# Write each sequence to the output file
for (i in 1:nrow(utr_sequences)) {
  # Write the gene name
  writeLines(paste(">", utr_sequences$ID[i], " | length =", utr_sequences$length[i], sep=""), utr_file_out)
  
  # Write the sequence
  writeLines(utr_sequences$FlippedSequence[i], utr_file_out)
}

# Close the output file
close(utr_file_out)



# entire genes (end to end) from 5'UTR to 3'UTR
# Define the output fasta file path
gene_output_file <- "Pf_genes.fasta"

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




#Plotting
#length_dist <- hist(utr_sequences$length)
library(ggplot2)
library(dplyr)

#min <- 0
#max <- 4800
ggplot(utr_sequences, aes(length)) + 
  geom_histogram(bins = 50, color="black", fill="grey") +
  scale_x_continuous(breaks = seq(0, max(utr_sequences$length), by = 200)) +
  geom_vline(xintercept = 1000, linetype = "dashed", color = "black") +  
  labs(title="Distribution of Pf annotated 5' UTR lengths (bp)") +
  theme_classic()



#Plotting motif frequencies
#AP2 motifs
genes_with_AP2_motifs <- read.csv("~/genes_with_AP2_motifs.csv", header=TRUE)


# categories reordered based on desc order of value
genes_with_AP2_motifs$y_label <- paste(genes_with_AP2_motifs$Domain, genes_with_AP2_motifs$AP2_binding_site, sep = " | ")
genes_with_AP2_motifs$y_label <- reorder(genes_with_AP2_motifs$y_label, genes_with_AP2_motifs$no_of_genes)
 
ggplot(genes_with_AP2_motifs, aes(x = no_of_genes, y = y_label)) +
  geom_col() + labs(title="AP2 domain core binding site", 
                    x = "Number of Genes", 
                    y = "ApiAP2") +
  theme_classic()


# # Find the maximum value, excluding the highest value
# max_value <- max(genes_with_AP2_motifs$no_of_genes)
# 
# # Filter the dataframe to exclude the row with the highest value
# genes_with_AP2_motifs_filtered <- genes_with_AP2_motifs[genes_with_AP2_motifs$no_of_genes != max_value, ]
# 
# 
# ggplot(genes_with_AP2_motifs_filtered, aes(y=AP2_binding_site,
#                                   x=no_of_genes)) + 
#   geom_col() + theme_classic()


#rsat motifs
genes_with_identified_motifs <- read.csv("~/genes_with_identified_motifs.csv", header=TRUE)


# Create a new factor variable with the categories reordered based on value
genes_with_identified_motifs$y_label <- paste(genes_with_identified_motifs$Motif, genes_with_identified_motifs$core.sequence, sep = " | ")
genes_with_identified_motifs$y_label <- reorder(genes_with_identified_motifs$y_label, genes_with_identified_motifs$no_of_genes)

ggplot(genes_with_identified_motifs, aes(y=y_label,
                                  x=no_of_genes)) + 
  geom_col() + labs(title="de novo motifs core binding site", 
                    x = "Number of Genes", 
                    y = "De novo motifs") + theme_classic()





# # Find the maximum value, excluding the highest value
# max_value2 <- max(genes_with_identified_motifs$no_of_genes)

# # Filter the dataframe to exclude the row with the highest value
# genes_with_identified_motifs_filtered <- genes_with_identified_motifs[genes_with_identified_motifs$no_of_genes != max_value2, ]
# 
# 
# ggplot(genes_with_identified_motifs_filtered, aes(y=identified_motifs,
#                                            x=no_of_genes)) + 
#   geom_col() + theme_classic()


#no stat sig diff
sig_test <- t.test(genes_with_AP2_motifs$no_of_genes, genes_with_identified_motifs$no_of_genes)


# repeat plotting for pwms
genes_with_AP2_pwms <- read.csv("~/genes_with_AP2_pwms.csv", header=TRUE)


# Create a new factor variable with the categories reordered based on value
genes_with_AP2_pwms$motif <- reorder(genes_with_AP2_pwms$motif, 
                                                  genes_with_AP2_pwms$no_of_matches)

ggplot(genes_with_AP2_pwms, aes(y=motif,
                                  x=no_of_matches)) + 
  geom_col() + labs(title="AP2 pwms") + theme_classic()





# # Find the maximum value, excluding the highest value
# max_value <- max(genes_with_AP2_pwms$no_of_genes)
# 
# # Filter the dataframe to exclude the row with the highest value
# genes_with_AP2_pwms_filtered <- genes_with_AP2_pwms[genes_with_AP2_pwms$no_of_genes != max_value, ]
# 
# 
# ggplot(genes_with_AP2_pwms_filtered, aes(y=AP2_binding_site,
#                                            x=no_of_genes)) + 
#   geom_col() + theme_classic()



genes_with_identified_pwms <- read.csv("~/genes_with_identified_pwms.csv", header=TRUE)


# Create a new factor variable with the categories reordered based on value
genes_with_identified_pwms$motif <- reorder(genes_with_identified_pwms$motif, 
                                                          genes_with_identified_pwms$no_of_matches)

ggplot(genes_with_identified_pwms, aes(y=motif,
                                         x=no_of_matches)) + 
  geom_col() + labs(title="de novo motifs pwms") + theme_classic()


#stat sig diff. may be due to longer length of AP2 pwms
sig_test_pwms <- t.test(genes_with_AP2_pwms$no_of_matches, genes_with_identified_pwms$no_of_matches)

# # Find the maximum value, excluding the highest value
# max_value2 <- max(genes_with_identified_pwms$no_of_genes)
# 
# # Filter the dataframe to exclude the row with the highest value
# genes_with_identified_pwms_filtered <- genes_with_identified_pwms[genes_with_identified_pwms$no_of_genes != max_value2, ]
# 
# 
# ggplot(genes_with_identified_pwms_filtered, aes(y=identified_pwms,
#                                                   x=no_of_genes)) + 
#   geom_col() + theme_classic()
