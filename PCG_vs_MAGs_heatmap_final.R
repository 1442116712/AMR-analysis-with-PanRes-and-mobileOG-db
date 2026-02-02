# install ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
# install.packages("vegan")
library(ComplexHeatmap)
library(circlize)
library(vegan)

# Nucleotide ####
# input data
data2 <- read.table('C:/Users/KEN/Desktop/UK/PhD/AMR/Rumen_previous/Combine_Panres_Heatmap/Pure_MAGs_updated/MAGs_Class.txt',
                   header = TRUE,sep = "\t", row.names = 1)
# change the data frame into array matrix
data_matrix2 <- as.matrix(data2) 

# Normalisation
dt2 <- decostand(t(data_matrix2), method = "total")
df2 <- t(dt2)

# Side annotation
rowsum2 <- read.table('C:/Users/KEN/Desktop/UK/PhD/AMR/Rumen_previous/Combine_Panres_Heatmap/Pure_MAGs_updated/MAGs_Class_sum.txt',
                      header = TRUE,sep = "\t", row.names = 1)

bar4 = rowAnnotation(
  ARG_abundance = anno_barplot(rowsum2,
                                  baseline = 0, bar_width = 0.9, width = unit(1.3, "cm"), gp = gpar(col = "white", fill = "#20B2AA"), 
                                  border = F, border_gp = gpar(lwd = 2), 
                                  axis_param = list(side = "bottom", at = c(0,150,300), labels = c("0","150","300"))), # direction = "reverse" - can reverse the direction 
  show_annotation_name = T, annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "bottom", annotation_name_rot = 90)


# Bottom annotation
phylum_colors2 <- c(rep("#FF69B4",2), rep("#FF8C00",1),rep("#000080",1),rep("#7FFFAA",3),rep("#FFD700",2),rep("#9932CC",2))
dummy_data2 <- rep(1, length(phylum_colors2)) 
ARG_type2 <- read.table('C:/Users/KEN/Desktop/UK/PhD/AMR/Rumen_previous/Combine_Panres_Heatmap/Pure_MAGs_updated/MAGs_Class_type.txt',header = TRUE,sep = "\t", row.names = 1)

bar5 <- HeatmapAnnotation(
  show_annotation_name = F,annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "right",annotation_name_rot = 0,
  . = anno_barplot(
    x = dummy_data2,
    gp = gpar(fill = phylum_colors2,col = NA), border = F,
    bar_width = 1, height = unit(0.2, "cm"),axis_param = list(at = "")),
  ARG_diversity = anno_barplot(    
    ARG_type2, 
    baseline = 0, bar_width = 0.9, height = unit(1.5, "cm"), gp = gpar(col = "white", fill = "#6495ED"), # bar setting
    border = F, border_gp = gpar(lwd = 2), 
    axis_param = list(direction = "reverse", side = "left", at = c(0,7,14), # axis parameter direction = "reverse", 
                      labels = c("")),ylim = c(0, 14)), gap = unit(1.5, "mm"))


# Top annotation
prevalence2 <- read.table('C:/Users/KEN/Desktop/UK/PhD/AMR/Rumen_previous/Combine_Panres_Heatmap/Pure_MAGs_updated/MAGs_Class_bar.txt',header = TRUE,sep = "\t", row.names = 1)
prevalence2$Percentage <- (prevalence2$ARG_detect / prevalence2$Metadata) * 100
prevalence_matrix2 <- as.matrix(prevalence2) 
Metadata2 <- prevalence_matrix2[ ,1]
ARG_detect2 <- prevalence_matrix2[ ,2]

# Row annotation
bar6 <- HeatmapAnnotation(
  ARG_prevalence = anno_barplot(
    cbind(prevalence2$Metadata, prevalence2$ARG_detect),
    gp = gpar(fill = c("#3CB371", "#FF69B4")), 
    beside = TRUE,                        
    bar_width = 0.85, height = unit(1.8, "cm"),
    border = F,
  ),
  show_annotation_name = F, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_name_side = "left", 
  annotation_name_rot = 0,
  
  Percentage = anno_text(
    paste0(round(prevalence2$Percentage, 0), "%"), 
    location = 0.5,  
    just = "center",            
    rot = 0,                   
    gp = gpar(fontsize = 6, fontface = "bold")     
  )
)

# Heatmap split
col_split2 <- c(rep("Culturable", 9), rep("Unculturable", 2))

col_fun = circlize::colorRamp2(c(0, 0.5, 1), c("#000080", "white", "red")) # the range should match the number of the colour

# heatmap setting
Heatmap2 = Heatmap(df2,
                   col = col_fun,
                   name = "Class in MAGs Nucleotide",
                   width = unit(6.1, "cm"), 
                   height = unit(8, "cm"), # size of the well # heatmap_width (including the text and title)
                   cluster_columns = F, 
                   cluster_rows = F, # data cluster
                   show_row_dend = F,
                   show_column_dend = F, 
                   row_dend_side = "left", # need to cluster first
                   column_dend_height = unit(0.5, "cm"), 
                   row_dend_width = unit(0.5, "cm"), # dendrogram height and location
                   column_title = "MAGs Nucleotide", 
                   row_title = "AMR Categories",
                   row_title_side = "right", 
                   column_title_side = "bottom", 
                   column_title_rot = F, # title location and rotation
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                   row_title_gp = gpar(fontsize = 10, fontface = "bold"), # fontface = "bold", fill = "green", col = "black", border = "black"), # title setting
                   right_annotation =bar4, 
                   bottom_annotation = bar5,
                   top_annotation = bar6,  # annotation chart
                   row_names_side = "right", 
                   column_names_side = "bottom", 
                   row_names_gp = gpar(fontsize = 8, fontface = "bold"), 
                   column_names_gp = gpar(fontsize = 8, fontface = "bold"), # location and size of the label name
                   show_heatmap_legend = F, 
                   show_column_names = T, 
                   show_row_names = F, 
                   column_split = col_split2, 
                   column_gap = unit(0.7, 'mm')) 
Heatmap2


# Heatmap3
# input data
data3 <- read.table('C:/Users/KEN/Desktop/UK/PhD/AMR/Rumen_previous/Combine_Panres_Heatmap/Pure_MAGs_updated/Pure_Class.txt',
                    header = TRUE,sep = "\t", row.names = 1)
# change the data frame into array matrix
data_matrix3 <- as.matrix(data3) 

# Normalisation
dt3 <- decostand(t(data_matrix3), method = "total")
df3 <- t(dt3)

# Side annotation
rowsum3 <- read.table('C:/Users/KEN/Desktop/UK/PhD/AMR/Rumen_previous/Combine_Panres_Heatmap/Pure_MAGs_updated/Pure_Class_sum.txt',
                      header = TRUE,sep = "\t", row.names = 1)

# bar chart annotation sum by column and row
bar7 = rowAnnotation(
  ARG_abundance = anno_barplot(rowsum3,
                             baseline = 0, bar_width = 0.9, width = unit(1.5, "cm"), gp = gpar(col = "white", fill = "#20B2AA"), 
                             border = F, border_gp = gpar(lwd = 2), 
                             axis_param = list(at = c(0,150,300), side = "bottom",labels = c("0","150","300"))), # direction = "reverse" - can reverse the direction 
  show_annotation_name = F, annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "bottom", annotation_name_rot = 90)


# Bottom annotation
phylum_colors3 <- c(rep("#FF69B4",2), rep("#FF8C00",1),rep("#DDA0DD",1),rep("#800000",1),rep("#000080",1),rep("#7FFFAA",4),rep("#FFD700",3))
dummy_data3 <- rep(1, length(phylum_colors3)) 
ARG_type3 <- read.table('C:/Users/KEN/Desktop/UK/PhD/AMR/Rumen_previous/Combine_Panres_Heatmap/Pure_MAGs_updated/Pure_Class_type.txt',header = TRUE,sep = "\t", row.names = 1)

bar8 <- HeatmapAnnotation(
  show_annotation_name = T,annotation_name_gp = gpar(fontsize = 8), annotation_name_side = "left",annotation_name_rot = 0,
  . = anno_barplot(
    x = dummy_data3,
    gp = gpar(fill = phylum_colors3,col = NA), border = F,
    bar_width = 1, height = unit(0.2, "cm"),axis_param = list(at = "")),
  ARG_diversity = anno_barplot(    
    ARG_type3, 
    baseline = 0, bar_width = 0.9, height = unit(1.5, "cm"), gp = gpar(col = "white", fill = "#6495ED"), # bar setting
    border = F, border_gp = gpar(lwd = 2), 
    axis_param = list(direction = "reverse", side = "left", at = c(0,7,14), # axis parameter direction = "reverse", 
                      labels = c("0","7","14")),ylim = c(0, 14)), gap = unit(1.5, "mm"))


# Top annotation
prevalence3 <- read.table('C:/Users/KEN/Desktop/UK/PhD/AMR/Rumen_previous/Combine_Panres_Heatmap/Pure_MAGs_updated/Pure_Class_bar.txt',header = TRUE,sep = "\t", row.names = 1)
prevalence3$Percentage <- (prevalence3$ARG_detect / prevalence3$Metadata) * 100
prevalence_matrix3 <- as.matrix(prevalence3) 
Metadata3 <- prevalence_matrix3[ ,1]
ARG_detect3 <- prevalence_matrix3[ ,2]

bar9 <- HeatmapAnnotation(
  ARG_prevalence = anno_barplot(
    cbind(prevalence3$Metadata, prevalence3$ARG_detect),
    gp = gpar(fill = c("#3CB371", "#FF69B4")), 
    beside = TRUE,                        
    bar_width = 0.85, height = unit(1.8, "cm"),
    border = F,
  ),
  show_annotation_name = T, 
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_name_side = "left", 
  annotation_name_rot = 0,

  Percentage = anno_text(
    paste0(round(prevalence3$Percentage, 0), "%"), 
    location = 0.5,  
    just = "center",            
    rot = 0,                   
    gp = gpar(fontsize = 6, fontface = "bold")     
  )
)

# heatmap setting
Heatmap3 = Heatmap(df3,
                   col = col_fun,
                   name = "Pure_class",
                   width = unit(7, "cm"), 
                   height = unit(8, "cm"), # size of the well # heatmap_width (including the text and title)
                   cluster_columns = F, 
                   cluster_rows = F, # data cluster
                   show_row_dend = F, 
                   show_column_dend = F, 
                   row_dend_side = "left", # need to cluster first
                   column_dend_height = unit(0.5, "cm"), 
                   row_dend_width = unit(0.5, "cm"), # dendrogram height and location
                   column_title = "PCGs Nucleotide", 
                   row_title = "AMR Categories",
                   row_title_side = "left", 
                   column_title_side = "bottom", 
                   column_title_rot = F, # title location and rotation
                   column_title_gp = gpar(fontsize = 10, fontface = "bold"), 
                   row_title_gp = gpar(fontsize = 10, fontface = "bold"), # fontface = "bold", fill = "green", col = "black", border = "black"), # title setting
                   right_annotation =bar7, 
                   bottom_annotation = bar8, 
                   top_annotation = bar9, # annotation chart
                   row_names_side = "left", 
                   column_names_side = "bottom", 
                   row_names_gp = gpar(fontsize = 8, fontface = "bold"), 
                   column_names_gp = gpar(fontsize = 8, fontface = "bold"), # location and size of the label name
                   show_heatmap_legend = F, 
                   show_column_names = T, 
                   show_row_names = T) 

Heatmap3


Heatmap0 = Heatmap3 + Heatmap2

Heatmap0

# legend setting 
phylum_labels <- c("Actinomycetota","Bacteroidota","Campylobacterota","Desulfobacteriota","Fibrobacterota","Firmicutes","Pseudomonadota","Verrucomicrobiota")
label_colors <- c("#FF69B4", "#FF8C00","#DDA0DD","#800000","#000080", "#7FFFAA","#FFD700","#9932CC")

lgd1 = Legend(labels = phylum_labels[1:8], labels_gp = gpar(fontsize = 8, fontface = "bold", col = "black"),
              legend_gp = gpar(fill = label_colors[1:8],cex = 0.5, frot = 1.2), 
              title_gp = gpar(fontsize = 10, fontface = "bold", col = "black"),
              title = "Phylum",  title_position = "topcenter", border = T,
              gap = unit(1, "cm"), ncol = 2)

lgd2 = Legend(col_fun = col_fun, title = "Percentage distribution", title_position = "topcenter", grid_width = unit(0.4, "cm"),legend_height = unit(1.9, "cm"),
              title_gp = gpar(fontsize = 10, fontface = "bold", col = "black"), border = "darkgrey")

#lgd3 = Legend(col_fun = col_fun1, title = "Percentage in pointchart", title_position = "topcenter", grid_width = unit(0.4, "cm"),legend_height = unit(1.5, "cm"),
#title_gp = gpar(fontsize = 10, fontface = "bold", col = "black"), border = "darkgrey")
pd = packLegend(lgd2, lgd1,column_gap = unit(3, "cm"),max_height = unit(25, "cm"), direction = "horizontal")
Heatmap0
draw(Heatmap0, ht_gap = unit(5, "mm"))
draw(pd, x = unit(20, "cm"), y = unit(26.2, "cm")) # 12 inched x 15 inches
draw(pd, x = unit(12, "cm"), y = unit(22, "cm")) # 850 Width x 900 height
