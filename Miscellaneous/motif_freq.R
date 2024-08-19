install.packages('seqinr') # do this once
library(seqinr)
library(ggplot2)
library(gridExtra)

rm(list=ls()) #remove objects from specified envt  

processFile = function(multifasta) {
  genes <- c()
  seqs <- c()
  con = file(multifasta, "r")
  sequence=""
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    print(line)
    ## header
    if(length(grep('>', line))){
      g=unlist(strsplit(line, split=" "))[1]
      g= gsub('>','',g) ### Remove '>'
      genes=c(genes,g)
      seqs=c(seqs, sequence)
      sequence=""
    }else{
      ### sequence
      sequence=paste0(sequence,line)
    }
    
  }
  seqs=c(seqs, sequence) ## last sequence
  #close(con)
  seqs=seqs[2:length(seqs)]
  dat=cbind(genes,seqs)
  return(dat)
}


# Create a list of cluster files
clusters_list <- list("allupreg8_cluster1pfpvp_filtered_251021.txt", "allupreg8_cluster2pfpvp_filtered_251021.txt", "allupreg8_cluster3pfpvp_filtered_251021.txt")

# Apply the processFile function to the clusters_list using lapply()
clusters_data <- lapply(clusters_list, processFile)




motif_list <-list(cluster1= c('AGACAA','TGTTGTA','AATTTTTAAA','ATTACA','CAACCAA','TGCATGCA'),
                  cluster2= c('TTATATAATA','AGGTAA','GCACTA'), 
                  cluster3=c('ATAGACAAA','ATGCACACA','AGCTAGCTA'))

for (cluster_name in names(motif_list)) {
  cluster_motifs <- motif_list[[cluster_name]]
  for (i in 1:length(cluster_motifs)){
    assign(paste0("frequencyMotif_", cluster_name, "_", toString(i)), 0) ## create a new variable
    motif <- tolower(cluster_motifs[i])
    l.motif <- nchar(motif)
    frequency <- c()
    for (k in seq_along(clusters_data)) {
      for(j in 1:dim(clusters_data[[k]])[1]){
        count <- 0
        gene_seq <- tolower(clusters_data[[k]][j,2])
        l.gene_seq <- nchar(gene_seq)
        for(ii in 1:((l.gene_seq-l.motif)+1)){
          test <- substr(gene_seq, ii, (ii+l.motif-1))
          cat(test , " .... ", nchar(test),"\n")
          cat(motif , " .... ", nchar(motif),"\n")
          if(tolower(motif) == tolower(test)){
            count <- count + 1
          }
        }
        frequency <- c(frequency, count)
      }
    }
    assign(paste0("frequencyMotif_", cluster_name, "_", toString(i)), frequency)
  }
}


freq_list <- list(
  frequencyMotif_1 = frequencyMotif_1,
  frequencyMotif_2 = frequencyMotif_2,
  frequencyMotif_3 = frequencyMotif_3,
  frequencyMotif_4 = frequencyMotif_4,
  frequencyMotif_5 = frequencyMotif_5,
  frequencyMotif_6 = frequencyMotif_6,
  frequencyMotif_7 = frequencyMotif_7,
  frequencyMotif_8 = frequencyMotif_8,
  frequencyMotif_9 = frequencyMotif_9,
  frequencyMotif_10 = frequencyMotif_10,
  frequencyMotif_11 = frequencyMotif_11,
  frequencyMotif_12 = frequencyMotif_12
)


freqresultlist <- list()
for (f in 1:length(freq_list)) {
  freq_matrix = as.matrix(freq_list[[f]])
  colnames(freq_matrix) = "freq"
  motif_freq = as.data.frame(cbind(dat, freq_matrix))
  # name = paste("dat",freq_list[f],sep="")
  freqresultlist[[f]] <- motif_freq
  
} 

name <- as.character(names(freq_list))
names(freqresultlist) <- name
