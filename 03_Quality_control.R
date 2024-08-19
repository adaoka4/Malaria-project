library(DESeq2) 
library(pheatmap)
library(viridis)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)


files <- file.path("~/test/", "PbPfPv_counts.csv")
PbPfPv_counts <- read.csv(files, header = T, sep = ",")
files2 <- file.path("~/count_data/", "metadata_PbPfPv_260421.csv")
Metadata_pbpfpv <- read.csv(files2, header = T, sep = ",")
Metadata_pbpfpv <- Metadata_pbpfpv %>% mutate_all(factor)


Countdata<- data.matrix(PbPfPv_counts[,-1])
rownames(Countdata) <- PbPfPv_counts$Pb_gene_stable_ID

dds <- DESeqDataSetFromMatrix(countData = Countdata, colData = Metadata_pbpfpv, design = ~parasite, ignoreRank = T)
nrow(dds)

vst <- vst(dds)
assay_vst <- assay(vst)

#pca on normal data (not batch corrected)
#with assay_vst
ntop <- 500
rv <- rowVars(assay_vst)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay_vst[select, ] )
# #mat2 <- limmavst[select, ]
pca2<-prcomp(mat)
View(pca2$x)
pca2_df <- as.data.frame(pca2$x)
View(pca2_df)


#determine % variance per PC
library("FactoMineR")
res.pca <- prcomp(mat)
library("factoextra")
eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_pca_ind(res.pca, col.ind = "black")


#plot with ggplot
pcaData <- data.frame(xPC=pca2_df$PC1, 
                      yPC=pca2_df$PC2, 
                      Parasite = Metadata_pbpfpv$parasite,
                      Species.Study = Metadata_pbpfpv$batch,
                      name= as.character(Metadata_pbpfpv$id),
                      stringsAsFactors = F)

pcaplot<- ggplot(pcaData, aes(x=xPC, y=yPC, colour = Parasite, shape = Species.Study, label=name)) 
# 
# 
pcaplot + geom_point() + scale_shape_manual(values = c(1:13)) +
#   #geom_text_repel(size = 2, force = 20, max.overlaps = Inf) + 
  xlab("PC1: 37.1% variance") +
 ylab("PC2: 20.9% variance") +# ggtitle("Before Batch effect Removal")+
 theme_bw() +
theme(panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size=0.5))
# 

######################################
#heatmap of dataset correlation on normal dataset (not batch corrected)
VarGenes <- rowVars(assay_vst) #batch corrected
#summary(VarGenes)
#q90 <- quantile(VarGenes, .90)
#limmavst_90 <- limmavst[VarGenes > q90,]
#dim(limmavst_90)
#Cor_limmavst_90 <- cor(limmavst_90, method = "spearman")

#pheatmap(Cor_limmavst_90, color = inferno(10), border_color = NA, main = "90th percentile",
#cluster_rows=TRUE, show_rownames=T, cluster_cols=TRUE, fontsize = 6)


top500 <- head(order(VarGenes, decreasing=TRUE), 500)
Count_top500 <- assay_vst[top500,]
cor_top500 <- cor(Count_top500, method = "spearman")

##pheatmap(cor_top500, color = inferno(7), 
#border_color = NA, main = "COR",
#cluster_rows=T, show_rownames=T, show_colnames = F, 
#cluster_cols=T, fontsize = 6, fontsize_row = 5)



#cols <- viridis(10,1,0,1,1,"inferno")
cols2 <- colorRamp2(breaks = c(-0.2,0,0.2,0.4,0.6,0.8,1),
                    colors = c("#1B0C42FF", "#4B0C6BFF", "#781C6DFF", "#A52C60FF",
                               "#CF4446FF", "#FB9A06FF", "#FCFFA4FF"), transparency = 0, space = "LAB")


heat <- Heatmap(cor_top500, col = cols2, name = "Dataset Correlation",
                color_space = "LAB", border = NA,
                cluster_rows = T, cluster_columns = T, 
                row_names_gp = gpar(fontsize = 4),
                show_column_names = F, show_row_names = T, 
                show_column_dend = F,
                show_heatmap_legend = T,
                heatmap_legend_param = list(
                  at = c(-0.2,0,0.2,0.4,0.6,0.8,1),
                  color_bar = 'continuous',
                  legend_direction = 'horizontal',
                  legend_width = unit(5, 'cm'),
                  legend_height = unit(2, 'cm'),
                  title_position = 'topcenter',
                  title_gp=gpar(fontsize = 8, fontface = 'bold'),
                  labels_gp=gpar(fontsize = 6)))

png("spearmancorplot.png", width = 4, height = 6, units = 'in', res = 500)
draw(heat, heatmap_legend_side = 'top')
dev.off()


############################
#batch effect correction
#with combat-seq
library(sva)
groups = sapply(as.character(Metadata_pbpfpv$parasite), switch, "Sporozoite" = 1, "Liver stages" = 2, "Asexual blood stages" = 3, "Gametocyte" = 4, "Ookinete" = 5, USE.NAMES = F)
batches = sapply(as.character(Metadata_pbpfpv$batch), switch, "Pb.Caldelari" = 1, "Pb.ToroMoreno" = 2, "Pb.Otto" = 3, "Pb.Yeoh" = 4, "Pb.Pandey" = 5, "Pb.Modryznska" = 6, "Pf.Zanghi" = 7,
                 "Pf.Lindner" = 8, "Pf.Shaw" = 9, "Pf.Fan" = 10, "Pf.Unpublished" = 11, "Pv.Muller" = 12, "Pv.Zhu"= 13, USE.NAMES = F)
combatcount <- ComBat_seq(counts = as.matrix(Countdata), batch = batches, group = groups) 
#combatvst <- ComBat_seq(counts = as.matrix(limmavst), batch = batches, group = groups) #gives Error: Negative counts not allowed, limma takes neg values

#running DESeq2 with combatcount as input from beginning before normalization and vst = PCA_combatcount_input
combatdds <- DESeqDataSetFromMatrix(countData = combatcount, colData = Metadata_pbpfpv, design = ~parasite, ignoreRank = T)
nrow(combatdds)


#dds <- DESeqDataSetFromMatrix(countData = Countdata, colData = Metadata_pbpfpv, design = ~batch+stage, ignoreRank = T)
combatdds <- combatdds[rowSums(counts(combatdds)) > 1, ]
nrow(combatdds)

#combatdds <- DESeq(combatdds) #not needed here

#plotDispEsts(combatdds)

combatvst <- vst(combatdds)

# plotPCA(vst, intgroup = "stage")
# plotPCA(vst, intgroup = "parasite")

Combat_vst <- assay(combatvst)

#remove batch effect with limma
#library(limma)
#limmavst <- removeBatchEffect(limmavst, batch = Metadata_pbpfpv$batch, design=model.matrix(~Metadata_pbpfpv$parasite))
#limmacount <- removeBatchEffect(Countdata, batch = Metadata_pbpfpv$batch, design=model.matrix(~Metadata_pbpfpv$parasite))


ntop <- 500
rv <- rowVars(Combat_vst)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( Combat_vst[select, ] )
# mat2 <- combatcount[select, ]
pca2<-prcomp(mat)
View(pca2$x)
pca2_df <- as.data.frame(pca2$x)
View(pca2_df)


#determine % variance per PC
library("FactoMineR")
res.pca <- prcomp(mat)
library("factoextra")
eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_pca_ind(res.pca, col.ind = "black")


#plot with ggplot
pcaData <- data.frame(xPC=pca2_df$PC1, 
                      yPC=pca2_df$PC2, 
                      Parasite = Metadata_pbpfpv$parasite,
                      Species.Study = Metadata_pbpfpv$batch,
                      name= as.character(Metadata_pbpfpv$id),
                      stringsAsFactors = F)

pcaplot<- ggplot(pcaData, aes(x=xPC, y=yPC, colour = Parasite, shape = Species.Study, label=name)) 

png("pcaplotnoBE.png", width = 6, height = 6, units = 'in', res = 500)
pcaplot + geom_point() + scale_shape_manual(values = c(1:13)) +
  #geom_text_repel(size = 2, force = 20, max.overlaps = Inf) + 
  xlab("PC1: 54% variance") +
  ylab("PC2: 21% variance") +# ggtitle("After Batch effect Removal")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))



#draw(heat, heatmap_legend_side = 'top')
dev.off()



#batch effect correction
#with limma
#remove batch effect with limma
library(limma)
limmavst <- removeBatchEffect(assay_vst, batch = Metadata_pbpfpv$batch, design=model.matrix(~Metadata_pbpfpv$parasite))
#limmacount <- removeBatchEffect(Countdata, batch = Metadata_pbpfpv$batch, design=model.matrix(~Metadata_pbpfpv$parasite))


ntop <- 500
rv <- rowVars(limmavst)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( limmavst[select, ] )
# mat2 <- combatcount[select, ]
pca2<-prcomp(mat)
View(pca2$x)
pca2_df <- as.data.frame(pca2$x)
View(pca2_df)


#determine % variance per PC
library("FactoMineR")
res.pca <- prcomp(mat)
library("factoextra")
eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_pca_ind(res.pca, col.ind = "black")


#plot with ggplot
pcaData <- data.frame(xPC=pca2_df$PC1, 
                      yPC=pca2_df$PC2, 
                      Parasite = Metadata_pbpfpv$parasite,
                      Species.Study = Metadata_pbpfpv$batch,
                      name= as.character(Metadata_pbpfpv$id),
                      stringsAsFactors = F)

pcaplot<- ggplot(pcaData, aes(x=xPC, y=yPC, colour = Parasite, shape = Species.Study, label=name)) 

png("pcaplotnoBE.png", width = 6, height = 6, units = 'in', res = 500)
pcaplot + geom_point() + scale_shape_manual(values = c(1:13)) +
  #geom_text_repel(size = 2, force = 20, max.overlaps = Inf) + 
  xlab("PC1: 62.2% variance") +
  ylab("PC2: 20.6% variance") +# ggtitle("After Batch effect Removal")+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.5))



#draw(heat, heatmap_legend_side = 'top')
dev.off()
