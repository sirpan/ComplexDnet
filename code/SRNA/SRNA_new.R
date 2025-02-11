library(Seurat)
library(clustree)

folder_list <- list.dirs(path = "D:/GSE182365_RAW", full.names = TRUE, recursive = FALSE)

seurat_list <- list()
samples <- list.files("D:/GSE182365_RAW", pattern = "_barcodes.tsv")
samples <- gsub("_barcodes.tsv.gz", "", samples)
data_dir <- "D:/GSE182365_RAW"

for (sample_id in samples) {
 
  target_dir <- file.path(data_dir, sample_id)
  if (!dir.exists(target_dir)) {
    dir.create(target_dir)
  }
  

  barcodes_file <- file.path(data_dir, paste0(sample_id, "_barcodes.tsv.gz"))
  features_file <- file.path(data_dir, paste0(sample_id, "_genes.tsv.gz"))  # 或 "_features.tsv"
  matrix_file <- file.path(data_dir, paste0(sample_id, "_matrix.mtx.gz"))
  
 
  file.rename(barcodes_file, file.path(target_dir, "barcodes.tsv.gz"))
  file.rename(features_file, file.path(target_dir, "features.tsv.gz"))
  file.rename(matrix_file, file.path(target_dir, "matrix.mtx.gz"))
}

for (folder in folder_list) {
 
  data <- Read10X(data.dir = folder)
  
 
  seurat_obj <- CreateSeuratObject(counts = data, project = basename(folder))
  
 
  seurat_list[[basename(folder)]] <- seurat_obj
}

data <- Read10X(data.dir = "D:/GSE182365_RAW/GSM5527640_HSD_hepatocytes")
seurat_obj <- CreateSeuratObject(counts = data, project = "GSM5527640_HSD_hepatocytes")
seurat_list["GSM5527640_HSD_hepatocytes"] <- seurat_obj
View(seurat_list)

seurat_merge <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
save(seurat_obj,file="D:\\metabolic network\\Wgcna\\SRNA_single_prepration.Rdata")
load("D:\\metabolic network\\Wgcna\\SRNA_single_prepration.Rdata")

library(dplyr)

library(Seurat)

library(patchwork)


seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")###如果是小鼠要改成"^mt-"


VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot1 + plot2


seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seurat_obj <- NormalizeData(seurat_obj)

seurat_obj<- FindVariableFeatures(seurat_obj, selection.method = "vst")

top10 <- head(VariableFeatures(seurat_obj), 10)

plot1 <- VariableFeaturePlot(seurat_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

# scale data
seurat_merge <- ScaleData(seurat_obj)

seurat_merge <- RunPCA(seurat_merge, features = VariableFeatures(object = seurat_merge))
save(seurat_merge,file="D:\\metabolic network\\Wgcna\\SRNA_PCA_single.Rdata")
load("D:\\metabolic network\\Wgcna\\SRNA_PCA_single.Rdata")

VizDimLoadings(seurat_merge, dims = 1:2, reduction = "pca")
DimPlot(seurat_merge, reduction = "pca") + NoLegend()
DimHeatmap(seurat_merge, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(seurat_merge)


seurat_merge <- FindNeighbors(seurat_merge, dims = 1:20)
my_colors <- c("#ea5c6f","#f7905a","#e187cb","#fb948d","#e2b159","#ebed6f","#b2db87","#7ee7bb","#64cccf","#a9dce6","#a48cbe","#e4b7d6"
)

seurat_merge <- FindClusters(seurat_merge, resolution = seq(0.1, 1.1, by = 0.1))
clustree_plot<-clustree(seurat_merge, prefix = "RNA_snn_res.",node_text_size = 5)
clustree_plot <- clustree_plot + scale_color_manual(values = my_colors) +
  theme(text = element_text(size = 20)  
  )



clustree_plot
ggsave(
  filename = "D:\\metabolic network\\Wgcna\\clustree_plot_new.png",  # 
  plot = clustree_plot,                                            # 
  width = 10,                                                      # 
  height = 8,                                                      # 
  dpi = 300                                                        # 
)

seurat_merge <- FindClusters(seurat_merge, resolution = 0.4)


seurat_merge <- RunUMAP(seurat_merge, dims = 1:20)
seurat_merge <- RunTSNE(seurat_merge, dims = 1:20)
colors<-c("#ea5c6f","#f7905a","#e187cb","#fb948d","#e2b159","#ebed6f","#b2db87","#7ee7bb","#64cccf","#a9dce6","#a48cbe","#e4b7d6")
colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C")

umap_plot <-DimPlot(seurat_merge, reduction = "umap",cols = colors)
umap_plot
ggsave(filename = "D:\\metabolic network\\Wgcna\\umap_plot.png", plot = umap_plot, width = 8, height = 6, dpi = 300)
DimPlot(seurat_merge, reduction = "tsne",cols = colors)
save(seurat_merge,file="D:\\metabolic network\\Wgcna\\SRNA_no_harmony.Rdata")

sample_colors <- c("#FF6F61", "#6B5B95", "#70C1B3", "#FFBE0B", "#FF85A1", "#A05195")
p <- DimPlot(seurat_merge, reduction = "tsne", group.by = "orig.ident",cols = sample_colors)
p
#harmony
seurat_merge_v5 <- IntegrateLayers(
  object = seurat_merge, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "cca",
  verbose = FALSE
)
#harmony
seurat_merge_v5 <- IntegrateLayers(
  object = seurat_merge, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

seurat_merge_v5 <- FindClusters(seurat_merge_v5, resolution = seq(0.1, 1.1, by = 0.1))



seurat_merge_v5 <- FindNeighbors(seurat_merge_v5, reduction = "cca", dims = 1:20)
seurat_merge_v5 <- FindClusters(seurat_merge_v5, resolution = 0.3, cluster.name = "cca_clusters")

seurat_merge_v5 <- RunUMAP(seurat_merge_v5, reduction = "cca", 
                           dims = 1:20, 
                           reduction.name = "umap.cca")
p5 <- DimPlot(
  seurat_merge_v5,
  reduction = "umap.cca",
  group.by = c("cca_clusters"),cols =colors)+
  DimPlot(
    seurat_merge_v5,
    reduction = "umap.cca",
    group.by = c("orig.ident"),cols =sample_colors)

p5



clustree(seurat_merge_v5, prefix = "RNA_snn_res.")

seurat_merge_v5 <- FindNeighbors(seurat_merge_v5, reduction = "harmony", dims = 1:20)

seurat_merge_v5 <- FindClusters(seurat_merge_v5, resolution = seq(0.1, 1.1, by = 0.1))
clustree(seurat_merge_v5, prefix = "RNA_snn_res.")
seurat_merge_v5 <- FindClusters(seurat_merge_v5, resolution = 0.2, cluster.name = "harmony_clusters")

seurat_merge_v5 <- RunUMAP(seurat_merge_v5, reduction = "harmony", 
                           dims = 1:20, 
                           reduction.name = "umap.harmony")
p7 <- DimPlot(
  seurat_merge_v5,
  reduction = "umap.harmony",
  group.by = c("harmony_clusters"),cols =colors)+
  DimPlot(
    seurat_merge_v5,
    reduction = "umap.harmony",
    group.by = c("orig.ident"),cols =sample_colors)

p7
save(seurat_merge_v5,file="D:\\metabolic network\\Wgcna\\SRNA_harmony_umap_0.2.Rdata")




#BiocManager::install("celldex")
#BiocManager::install("SingleR")
library(celldex)
library(SingleR)
library(BiocParallel)
library(ggplot2)
load(file="D:\\metabolic network\\Wgcna\\SRNA_no_harmony.Rdata")


Idents(seurat_merge) <- seurat_merge$cca_clusters
seurat_merge <- JoinLayers(seurat_merge, layers = c("counts", "data", "scale.data"))
allmarkers <- FindAllMarkers(seurat_merge,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
allmarkers <- FindAllMarkers(seurat_merge,only.pos = TRUE)
write.csv(allmarkers, file = "D:\\metabolic network\\Wgcna\\all_markers.csv", row.names = FALSE)

top5_markers <- allmarkers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
top10_markers <- allmarkers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

top5_markers


DotPlot(seurat_merge, features = unique(markergenes)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 8, hjust = 1)) + NoLegend()


DoHeatmap(seurat_merge, features = unique(top5_markers$gene)) +
  theme(axis.text.y = ggplot2::element_text(size = 8)) + NoLegend()


save(allmarkers, top5_markers, file = "markers_integrated.Rdata")


markergenes <- c("Ptprd",'Rorc','A2ml1','Hck','Cd274')

markergenes <- c("Cd274")
p3 <- FeaturePlot(seurat_merge, 
                  features = markergenes, 
                  cols = c("#CAD8DA", "#f57c6e")) +  # Set color gradient to red shades
  NoLegend()+
  theme(
    axis.title = element_blank(),      # Remove axis titles
    axis.text = element_blank(),       # Remove axis text labels
    axis.ticks = element_blank(),      # Remove axis ticks
    axis.line = element_blank(),       # Remove axis lines
    panel.grid = element_blank(),      # Remove grid lines (optional)
    panel.border = element_blank(),    # Remove panel borders (optional)
    panel.background = element_blank(),# Remove panel background (optional)
    plot.margin = unit(c(0, 0, 0, 0), "cm")  # Remove plot margins (optional)
  )

# Display the plot (optional)
print(p3)
ggsave(filename = "D:\\metabolic network\\Wgcna\\Cd274_plot.png", plot = p3, width = 6, height = 6, dpi = 300)

save(allmarkers, seurat_merge,seurat_obj,file = "markers_integrated.Rdata")


#dotplot图
library(Seurat)
library(ggplot2)

# Define a multi-color palette
custom_colors <- c('#f2b56f', '#f57c6e', '#84c3b7', '#b8aeed')

# Generate the DotPlot with customized size and multi-color gradient
dotplot_custom <- DotPlot(seurat_merge, 
                          features = unique(markergenes), 
                          dot.scale = 8,  # Adjust dot size
                          cols = custom_colors) +  # Initial color gradient
  scale_color_gradientn(colors = custom_colors) +  # Customize color bar with multiple colors
  scale_size(range = c(2, 12)) +  # Adjust dot size range
  theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1)) +  # Rotate x-axis labels
  guides(color = guide_colorbar(title = "Expression Level", 
                                barwidth = 1, 
                                barheight = 10)) +  # Customize color bar title and size
  theme(legend.position = "right",  # Position of the legend
        legend.title = element_text(size = 12, face = "bold"),  # Legend title style
        legend.text = element_text(size = 10))  # Legend text style

# Display the customized DotPlot
print(dotplot_custom)

# Save the plot if desired
ggsave(filename = "D:\\metabolic network\\Wgcna\\custom_dotplot_with_legend.png", plot = dotplot_custom, width = 10, height = 8, dpi = 300)
dotplot_custom <- DotPlot(seurat_merge, 
                          features = unique(markergenes), 
                          dot.scale = 8,  # Adjust dot size
                          cols = custom_colors) +  # Initial color gradient
  scale_color_gradientn(colors = custom_colors) +  # Customize color bar with multiple colors
  scale_size(range = c(2, 12)) +  # Adjust dot size range
  theme(axis.title.y = element_text(size = 16),  # Increase y-axis title size
        axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1)) +  # Rotate x-axis labels
  guides(color = guide_colorbar(title = "Expression Level", 
                                barwidth = 1, 
                                barheight = 10)) +  # Customize color bar title and size
  theme(legend.position = "right",  # Position of the legend
        legend.title = element_text(size = 12, face = "bold"),  # Legend title style
        legend.text = element_text(size = 10))  # Legend text style

# Display the customized DotPlot
print(dotplot_custom)

