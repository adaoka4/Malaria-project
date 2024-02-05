library(ggplot2)
library(reshape2)
library(readxl)
library(factoextra)
library(NbClust)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)

#load DEGs and make bar plot
DEG_summary <- read.csv("~/Pb_070421/DEG_summary_200521.csv", header = T, sep = ",")
DEG_summary <- melt(DEG_summary)
ggplot <- ggplot(DEG_summary, aes(x=sample, y=value, fill=variable))+ geom_bar(stat="identity", color="black")

png("DEGplot.png", width = 6, height = 6, units = 'in', res = 500)

ggplot + scale_x_discrete(limits=c("Sporozoite_Caldelari","Sporozoite_ToroMoreno",
                            "L6h_Caldelari","L24h_Caldelari","L24h_ToroMoreno",
                            "L36h_ToroMoreno","L48h_Caldelari","L48h_ToroMoreno",
                            "L54h_Caldelari","L60h_Caldelari",
                            "Ring_Otto","Trophozoite_Otto","Schizont_Otto","Asexual_BS_Yeoh",
                            "Gametocyte_Otto","Gametocyte_Pandey","Ookinete_Otto","Ookinete_Modrzynska"))+ 
labs(title="Differentially expressed genes per sample",x="Samples",y = "# DE Genes")+
scale_fill_manual(name = "", labels = c("Upregulated", "Downregulated"), values=c('#BD3039','#2E5894'))+
theme_classic()+
theme(axis.text.x=element_text(angle = 90, hjust = 1))

dev.off()


DEG_data_frame <- read_excel("~/Pb_070421/DEG_list_200521.xlsx")
DEGfinallist <- as.list(data_DEG_frame)

#make volcano plots
for (i in 1:length(DEGfinallist)) {
volcano <- DEGfinallist[[18]] 
#png(paste(i,"plot.png"), width = 6, height = 6, units = 'in', res = 500)
EnhancedVolcano(volcano,
                lab = rownames(volcano),
                x = 'log2FoldChange',
                y = 'padj', title = as.character(names(DEGfinallist))[18],
                subtitle = "", pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 1.0, col=c('grey', 'grey', 'grey', 'red3'),
                labSize = 2.0, caption = "", legendPosition = 'none')

#dev.off()

}

#transform data for quality control
vst <- vst(dds)
Count_vst <- assay(vst)

Countavrep <- cbind((Count_vst[,1]+Count_vst[,2]+Count_vst[,3]+Count_vst[,4]+Count_vst[,5]+Count_vst[,6])/6,
                    (Count_vst[,7]+Count_vst[,8])/2,
                    (Count_vst[,9]+Count_vst[,10]+Count_vst[,11]+Count_vst[,12]+Count_vst[,13]+Count_vst[,14]+Count_vst[,15]+Count_vst[,16]+Count_vst[,17]+Count_vst[,18])/10,
                    (Count_vst[,19]+Count_vst[,20]+Count_vst[,21]+Count_vst[,22])/4,
                    (Count_vst[,23]+Count_vst[,24]+Count_vst[,25]+Count_vst[,26]+Count_vst[,27]+Count_vst[,28])/6,
                    (Count_vst[,29]+Count_vst[,30])/2,
                    (Count_vst[,31]+Count_vst[,32])/2,
                    (Count_vst[,33]+Count_vst[,34]+Count_vst[,39]+Count_vst[,40]+Count_vst[,41])/5,
                    (Count_vst[,35]+Count_vst[,36]+Count_vst[,39]+Count_vst[,40]+Count_vst[,41])/5,
                    (Count_vst[,37]+Count_vst[,38]+Count_vst[,39]+Count_vst[,40]+Count_vst[,41])/5,
                    (Count_vst[,42]+Count_vst[,43]+Count_vst[,44]+Count_vst[,45])/4,
                    (Count_vst[,46]+Count_vst[,47]+Count_vst[,48]+Count_vst[,49]+Count_vst[,50])/5)

colnames(Countavrep) <- c("Sporozoite","L6h","L24h","L36h","L48h","L54h","L60h","Ring","Trophozoite","Schizont","Gametocyte","Ookinete")

allupreg8stages <- read.csv("~/Pb_070421/allupreg8stages.csv", header = F, sep = ",")
allupreg8stages <- as.matrix(allupreg8stages)
allupreg8stages_countsavrep <- Countavrep[allupreg8stages,]
allupreg8stages_countsavrep <- allupreg8stages_countsavrep[,-c(2,4,6:7)]
allupreg8stages_countsavrep_scaled <- t(scale(t(allupreg8stages_countsavrep)))
write.table(allupreg8stages_countsavrep_scaled, "allupreg8stagesclustcountsavrep.csv", sep = ",", col.names = NA)

#detemine optmal k value
#elbow plot of kmeans wrapper
fviz_nbclust(allupreg8stages_countsavrep_scaled, cluster::pam, method = "wss", nboot = 1000, nstart = 25) +
  labs(title = "All upreg8 optimal number of pam clusters")

#silhouette plot of kmeans wrapper
fviz_nbclust(allupreg8stages_countsavrep_scaled, cluster::pam, method = "silhouette", nboot = 1000, nstart = 25) +
  labs(title = "All upreg8 optimal number of pam clusters")

#nbclust
nbclustallupreg8 <- NbClust(data= allupreg8stages_countsavrep_scaled, 
                            distance = "euclidean", min.nc = 2,
                            max.nc = 7, method = "kmeans",
                            index=c("kl","ch","cindex","db_order_height","silhouette",
                                    "duda","pseudot2","beale","ratkowsky",
                                    "gap","gamma","mcclain","gplus",
                                    "tau","sdindex","sdb_order_heightw"))

#heatmap of gene expression
myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

pamClustersallupreg8 <- cluster::pam(allupreg8stages_countsavrep_scaled, k = 3) # pre-select k = 3 centers
pamClustersallupreg8$clustering <- paste0('Cluster ', pamClustersallupreg8$clustering)

# fix order of the clusters to have 1 to 3, top to bottom
pamClustersallupreg8$clustering <- factor(pamClustersallupreg8$clustering,
                                          levels = c('Cluster 1', 'Cluster 2', 'Cluster 3'))

hmapallupreg8 <- Heatmap(allupreg8stages_countsavrep_scaled,
                         
                         # split the genes / rows according to the PAM clusters
                         split = pamClustersallupreg8$clustering,
                         cluster_row_slices = FALSE,
                         
                         name = 'Gene\nZ-\nscore',
                         
                         col = colorRamp2(myBreaks, myCol),
                         
                         # parameters for the colour-bar that represents gradient of expression
                         heatmap_legend_param = list(
                           color_bar = 'continuous',
                           legend_direction = 'vertical',
                           legend_width = unit(8, 'cm'),
                           legend_height = unit(5.0, 'cm'),
                           title_position = 'topcenter',
                           title_gp=gpar(fontsize = 12, fontface = 'bold'),
                           labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                         
                         # row (gene) parameters
                         cluster_rows = TRUE,
                         show_row_dend = TRUE,
                         #row_title = 'Statistically significant genes',
                         row_title_side = 'left',
                         row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                         row_title_rot = 90,
                         show_row_names = FALSE,
                         row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                         row_names_side = 'left',
                         row_dend_width = unit(25,'mm'),
                         
                         # column (sample) parameters
                         cluster_columns = TRUE,
                         show_column_dend = TRUE,
                         column_title = '8 stages',
                         column_title_side = 'bottom',
                         column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                         column_title_rot = 0,
                         show_column_names = TRUE,
                         column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                         column_names_max_height = unit(10, 'cm'),
                         column_dend_height = unit(25,'mm'),
                         
                         # cluster methods for rows and columns
                         clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                         clustering_method_columns = 'ward.D2',
                         clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                         clustering_method_rows = 'ward.D2')


draw(hmapallupreg8,
     heatmap_legend_side = 'left',
     row_sub_title_side = 'left')   



#######cluster with sub-clusters
# Perform initial clustering
myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

pamClustersallupreg8 <- cluster::pam(allupreg8stages_countsavrep_scaled, k = 3) # pre-select k = 3 centers
pamClustersallupreg8$clustering <- paste0('Cluster ', pamClustersallupreg8$clustering)

# fix order of the clusters to have 1 to 3, top to bottom
pamClustersallupreg8$clustering <- factor(pamClustersallupreg8$clustering,
                                          levels = c('Cluster 1', 'Cluster 2', 'Cluster 3'))

# Specify the number of sub-clusters for each initial cluster
num_sub_clusters <- c(4, 5, 6)

# Initialize a list to store sub-cluster results
sub_cluster_labels <- list()

# Iterate through the initial clusters
for (i in 1:3) {
  # Extract genes in the current initial cluster
  cluster_genes <- allupreg8stages_countsavrep_scaled[pamClustersallupreg8$clustering == paste0('Cluster ', i), ]
  
  # Perform sub-clustering using k-means with the specified number of sub-clusters
  sub_cluster_labels[[i]] <- kmeans(cluster_genes, centers = num_sub_clusters[i])$cluster
  
  # Modify sub-cluster labels to distinguish from the initial clusters
  sub_cluster_labels[[i]] <- paste0('Cluster ', i, '.Sub', sub_cluster_labels[[i]])
}

# Combine the sub-cluster labels into a single vector
sub_cluster_labels <- unlist(sub_cluster_labels)

# Modify the labels of the original clustering to include sub-cluster labels
pamClustersallupreg8$clustering <- sub_cluster_labels

# Create the heatmap with modified clustering labels
hmapallupreg8 <- Heatmap(allupreg8stages_countsavrep_scaled,
                         
                         # split the genes / rows according to the PAM clusters
                         split = pamClustersallupreg8$clustering,
                         cluster_row_slices = FALSE,
                         
                         name = 'Gene\nZ-\nscore',
                         
                         col = colorRamp2(myBreaks, myCol),
                         
                         # parameters for the colour-bar that represents gradient of expression
                         heatmap_legend_param = list(
                           color_bar = 'continuous',
                           legend_direction = 'vertical',
                           legend_width = unit(8, 'cm'),
                           legend_height = unit(5.0, 'cm'),
                           title_position = 'topcenter',
                           title_gp=gpar(fontsize = 12, fontface = 'bold'),
                           labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                         
                         # row (gene) parameters
                         cluster_rows = TRUE,
                         show_row_dend = TRUE,
                         #row_title = 'Statistically significant genes',
                         row_title_side = 'left',
                         row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                         row_title_rot = 90,
                         show_row_names = FALSE,
                         row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                         row_names_side = 'left',
                         row_dend_width = unit(25,'mm'),
                         
                         # column (sample) parameters
                         cluster_columns = TRUE,
                         show_column_dend = TRUE,
                         column_title = '8 stages',
                         column_title_side = 'bottom',
                         column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                         column_title_rot = 0,
                         show_column_names = TRUE,
                         column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                         column_names_max_height = unit(10, 'cm'),
                         column_dend_height = unit(25,'mm'),
                         
                         # cluster methods for rows and columns
                         clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                         clustering_method_columns = 'ward.D2',
                         clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                         clustering_method_rows = 'ward.D2')

# Draw the heatmap
draw(hmapallupreg8,
     heatmap_legend_side = 'left',
     row_sub_title_side = 'left')
###########################################


#04-02-22
allupreg8stages_countsavrep_scaled <- read.csv("~/allupreg8stagesclustcountsavrep.csv", header = T, sep = ",")
cluster2_aggtaa <- read.csv("~/aggtaa_motif_locations_w_oldids.txt", header = T, sep = "\t")
allupreg8stages_countsavrep_scaled_as_df <- as.data.frame(allupreg8stages_countsavrep_scaled)
colnames(allupreg8stages_countsavrep_scaled_as_df)[1] <- "Ensembl_id"

extract_aggtaa <- semi_join(allupreg8stages_countsavrep_scaled_as_df,cluster2_aggtaa)
rownames(extract_aggtaa) <- extract_aggtaa[,1]
extract_aggtaa[,1] <- NULL

#clustering for aggtaa and gcacta motif target genes
hmapaggtaa <- Heatmap(extract_aggtaa,
                         
                         # split the genes / rows according to the PAM clusters
                         #split = pamClustersallupreg8$clustering,
                         #cluster_row_slices = FALSE,
                         
                         name = 'Z-score',
                         
                         col = colorRamp2(myBreaks, myCol),
                         
                         # parameters for the colour-bar that represents gradient of expression
                         heatmap_legend_param = list(
                           color_bar = 'continuous',
                           legend_direction = 'vertical',
                           legend_width = unit(8, 'cm'),
                           legend_height = unit(4.0, 'cm'),
                           title_position = 'topcenter',
                           title_gp=gpar(fontsize = 8, fontface = 'bold'),
                           labels_gp=gpar(fontsize = 8, fontface = 'bold')),
                         
                         # row (gene) parameters
                         cluster_rows = FALSE,
                         show_row_dend = FALSE,
                         #row_title = 'AGGTAA',
                         row_title_side = 'right',
                         row_title_gp = gpar(fontsize = 4,  fontface = 'bold'),
                         row_title_rot = 90,
                         show_row_names = TRUE,
                         row_names_gp = gpar(fontsize = 4, fontface = 'bold'),
                         row_names_side = 'right',
                         row_dend_width = unit(25,'mm'),
                         
                         # column (sample) parameters
                         cluster_columns = FALSE,
                         show_column_dend = FALSE,
                         #column_title = '8 stages',
                         column_title_side = 'bottom',
                         column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                         column_title_rot = 0,
                         show_column_names = TRUE,
                         column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                         column_names_max_height = unit(10, 'cm'),
                         column_dend_height = unit(25,'mm'),
                         
                         # cluster methods for rows and columns
                         #clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                         #clustering_method_columns = 'ward.D2',
                         #clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                         #clustering_method_rows = 'ward.D2')
                          )

draw(hmapaggtaa,
     heatmap_legend_side = 'left',
     row_sub_title_side = 'right')   

png("aggtaaheatplot.png", width = 8, height = 12, units = 'in', res = 1500)
draw(hmapaggtaa, heatmap_legend_side = 'left')
dev.off()


cluster2_gcacta <- read.csv("~/gcacta_motif_locations_w_oldids.txt", header = T, sep = "\t")
extract_gcacta <- semi_join(allupreg8stages_countsavrep_scaled_as_df,cluster2_gcacta)
rownames(extract_gcacta) <- extract_gcacta[,1]
extract_gcacta[,1] <- NULL

hclust_gcacta_stages <- hclust(dist(t(extract_gcacta)))
hclust_gcacta_genes <- hclust(dist(extract_gcacta))
plot(hclust_gcacta_genes)
plot(hclust_gcacta)

ordering_height <- c(4,5,6,2,3,7,8,1) #this ordering was determined manually from visual inspection of the dendrogram from ring (down) to sporozoite (up), not left to right which is 1,7,8,6,4,5,2,3
db_order_height=extract_gcacta[,ordering_height]
new_order_dist_genes <- extract_gcacta[test,]

mat_genes = hclust_gcacta_genes[["merge"]]
mat_stages = hclust_gcacta_stages[["merge"]]
#indexes=c(mat[,1], mat[,2])

dim(mat)
get_index = function(mat){
  indexes=c()
  for ( i in 1:dim(mat)[1]){
    for ( j in 1:2){
      item= mat[i,j]
      if(item<0){
        indexes=c(indexes, abs(item))
      }
      
    }
  }
  return(indexes)   
}

test = get_index(mat)
index_stages = get_index(mat_stages)
index_genes = get_index(mat_genes)

new_db= extract_gcacta[index_genes, ]
last_df=new_db[,index_stages]
new_df= dist_gcacta_genes[,test]

#gcacta
hmapdb_new_order_dist_genes <- Heatmap(new_order_dist_genes,
                      
                      # split the genes / rows according to the PAM clusters
                      #split = pamClustersallupreg8$clustering,
                      #cluster_row_slices = FALSE,
                      
                      name = 'Z-score',
                      
                      col = colorRamp2(myBreaks, myCol),
                      
                      # parameters for the colour-bar that represents gradient of expression
                      heatmap_legend_param = list(
                        color_bar = 'continuous',
                        legend_direction = 'vertical',
                        legend_width = unit(8, 'cm'),
                        legend_height = unit(4.0, 'cm'),
                        title_position = 'topcenter',
                        title_gp=gpar(fontsize = 8, fontface = 'bold'),
                        labels_gp=gpar(fontsize = 8, fontface = 'bold')),
                      
                      # row (gene) parameters
                      cluster_rows = FALSE,
                      show_row_dend = FALSE,
                      #row_title = 'GCACTA',
                      row_title_side = 'right',
                      row_title_gp = gpar(fontsize = 4,  fontface = 'bold'),
                      row_title_rot = 90,
                      show_row_names = TRUE,
                      row_names_gp = gpar(fontsize = 4, fontface = 'bold'),
                      row_names_side = 'right',
                      row_dend_width = unit(25,'mm'),
                      
                      # column (sample) parameters
                      cluster_columns = FALSE,
                      show_column_dend = FALSE,
                      #column_title = '8 stages',
                      column_title_side = 'bottom',
                      column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                      column_title_rot = 0,
                      show_column_names = TRUE,
                      column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                      column_names_max_height = unit(10, 'cm'),
                      column_dend_height = unit(25,'mm'),
                      
                      # cluster methods for rows and columns
                      #clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                      #clustering_method_columns = 'ward.D2',
                      #clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                      #clustering_method_rows = 'ward.D2')
)

draw(hmapdb_new_order_dist_genes,
     heatmap_legend_side = 'left',
     row_sub_title_side = 'right')

png("gcactaheatplot.png", width = 8, height = 12, units = 'in', res = 1500)
draw(hmapaggtaa, heatmap_legend_side = 'left')
dev.off()

dist_gcacta_genes <- as.matrix(dist(extract_gcacta))
dist_gcacta_stages <- as.matrix(dist(t(extract_gcacta)))


#add cluster assignment to data
allupreg8clustered <- data.frame(
  allupreg8 =rownames(allupreg8stages_countsavrep),
  cluster=pamClustersallupreg8$clustering)

write.table(allupreg8clustered, "allupreg8clustered.csv", sep = ",", col.names = NA)


allupreg8_cluster1 <- read.csv("~/allupreg8_cluster1.csv", header = F, sep = ",")
allupreg8_cluster1 <- as.matrix(allupreg8_cluster1)

library(genefilter)
allupreg8_cluster1_countsavrep_scaled <- allupreg8stages_countsavrep_scaled[allupreg8_cluster1,]
allupreg8_cluster1_varGenes <- rowVars(allupreg8_cluster1_countsavrep_scaled)

allupreg8_cluster1_varGenes2 <- head(order(rowVars(allupreg8_cluster1_countsavrep_scaled), decreasing = TRUE), 10) 
allupreg8_cluster1_top10vargenes <- allupreg8_cluster1_countsavrep_scaled[allupreg8_cluster1_varGenes2,]

q90 <- quantile(VarGenes, .90)
Count_rld_90 <- Count_rld[VarGenes < q90,]


#redo clusterinh with k=12
library(RColorBrewer)
myCol <- colorRampPalette(c('dodgerblue', 'black', 'yellow'))(100)
myBreaks <- seq(-3, 3, length.out = 100)

ha <- rowAnnotation(foo = anno_block(
                       gp = gpar(col = 'white', fill = my_color_palette, alpha = 0.7),
                       labels = c(10, 4, 11, 2, 9, 7, 8, 5, 6, 3, 12, 1), labels_rot = -90), 
                    empty = anno_empty(border = FALSE, height = unit(8, "mm")))


my_color_palette <- brewer.pal(12, "Paired")
pamClustersallupreg8 <- cluster::pam(allupreg8stages_countsavrep_scaled, k = 12) # pre-select k = 3 centers
#pamClustersallupreg8$clustering <- paste0('Cluster ', pamClustersallupreg8$clustering)

# fix order of the clusters to have 1 to 3, top to bottom
#pamClustersallupreg8$clustering <- factor(pamClustersallupreg8$clustering,
                                          #levels = c('Cluster 1', 'Cluster 2', 'Cluster 3'))

hmapallupreg8 <- Heatmap(allupreg8stages_countsavrep_scaled,
                         
                         # split the genes / rows according to the PAM clusters
                         split = pamClustersallupreg8$clustering,
                         row_gap = unit(0, "mm"),
                         cluster_row_slices = TRUE,
                         
                         name = 'Gene\nZ-\nscore',
                         
                         col = colorRamp2(myBreaks, myCol),
                         
                         # parameters for the colour-bar that represents gradient of expression
                         heatmap_legend_param = list(
                           color_bar = 'continuous',
                           legend_direction = 'vertical',
                           legend_width = unit(8, 'cm'),
                           legend_height = unit(5.0, 'cm'),
                           title_position = 'topcenter',
                           title_gp=gpar(fontsize = 12, fontface = 'bold'),
                           labels_gp=gpar(fontsize = 12, fontface = 'bold')),
                         
                         # row (gene) parameters
                         cluster_rows = TRUE,
                         #show_row_dend = TRUE,
                         row_title = NULL,
                         #row_title_side = 'left',
                         #row_title_gp = gpar(fontsize = 12,  fontface = 'bold'),
                         #row_title_rot = 90,
                         show_row_names = FALSE,
                         #row_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                         #row_names_side = 'left',
                         row_dend_width = unit(25,'mm'),
                         
                         # column (sample) parameters
                         cluster_columns = TRUE,
                         show_column_dend = TRUE,
                         column_title_side = 'bottom',
                         column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
                         column_title_rot = 0,
                         show_column_names = TRUE,
                         column_names_gp = gpar(fontsize = 10, fontface = 'bold'),
                         column_names_max_height = unit(10, 'cm'),
                         column_dend_height = unit(25,'mm'),
                         
                         # cluster methods for rows and columns
                         clustering_distance_columns = function(x) as.dist(1 - cor(t(x))),
                         clustering_method_columns = 'ward.D2',
                         clustering_distance_rows = function(x) as.dist(1 - cor(t(x))),
                         clustering_method_rows = 'ward.D2',
                         
                         right_annotation = ha)




draw(hmapallupreg8,
     heatmap_legend_side = 'left')
     #row_sub_title_side = 'left')  


#add annotations
library(GetoptLong)  # for the function qq()
group_block_anno = function(group, empty_anno, gp = gpar(), 
                            label = NULL, label_gp = gpar()) {
  
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gp)
  if(!is.null(label)) {
    grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp, rot = -90)
  }
}

group_block_anno(1:4, "empty", gp = gpar(col = 'white', fill = "blue", alpha = 0.5), label = "cluster 2")
group_block_anno(5:8, "empty", gp = gpar(col = 'white', fill = "red", alpha = 0.5), label = "cluster 3")
group_block_anno(9:12, "empty", gp = gpar(col = 'white', fill = "purple", alpha = 0.5), label = "cluster 1")

allupreg8_12clusters <- data.frame(
  allupreg8 =rownames(allupreg8stages_countsavrep),
  cluster=pamClustersallupreg8$clustering)

#write.table(allupreg8_12clusters, "allupreg8_12clusters.csv", sep = ",", col.names = NA)
