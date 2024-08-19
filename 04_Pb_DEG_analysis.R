library(DESeq2) 
library(dplyr)
library(purrr)
library(rlist)
library(EnhancedVolcano)

#Load count data and metadata
Pb_counts <- read.csv("~/count_data/Pb_counts_names.csv", header = T, sep = ",")
Metadata <- read.csv("~/count_data/metadata_Pb_070421.csv", header = T, sep = ",")
Metadata <- Metadata %>% mutate_all(factor)
#genelength <- read.csv("~/count_data/genelength.csv", header = T)

stage1 <- read.csv("C:/Users/user/Desktop/Pb_070421/stage.txt", header = F, sep = "/")
stage2<- stage1[order(stage1$V1, decreasing = T), ]
stage1 <- as.vector(stage1$V1)

#counts for differential analysis
Countdata<- data.matrix(Pb_counts[,-1])
rownames(Countdata) <- Pb_counts$Ensemble_Pb_gene_stable_ID

#Run deseq
dds <- DESeqDataSetFromMatrix(countData = Countdata, colData = Metadata, design = ~stage, ignoreRank = T)
dds <- dds[rowSums(counts(dds)) > 1, ]
#colSums(counts(dds))

#outlier filtering/replacement turned off because no outlier was replaced when turned on
dds <- DESeq(dds)#, minReplicatesForReplace=Inf)
lfc_threshold <- log2(4) #setting lfc threshold in result() gives a more statistically confident result
p_value_cutoff <- 0.01

#perform all possible pairwise comparisons https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html
resultlist <- list()

for (i in 1:length(stage1)) {
  for (j in 1:length(stage2)) {
    if (stage1[i] == stage2[j]) {
      next
    }
    result <- results(dds,
                      contrast = c("stage", stage1[i], as.character(stage2[j])),
                      alpha = p_value_cutoff,
                      lfcThreshold = lfc_threshold,
                      cooksCutoff=FALSE)
    name <- paste("res",stage1[i],"vs",stage2[j], sep='_')
    resultlist[[name]] <- result
  }
}

#extract DEGs per comparison
DEGlist <- map(1:length(resultlist), function(x) {
  y <- resultlist[[x]]
  resSig_lfc <- subset(y, padj < 0.01 & abs(log2FoldChange) > 2)
})

name <- as.character(names(resultlist))
names(DEGlist) <- name

#export summary stats for DEGs per comparison
for (x in 1:length(resultlist)) {
  sink(paste(as.character(names(resultlist))[x],".txt"))
  print(summary(resultlist[[x]]))
  sink()
 }
