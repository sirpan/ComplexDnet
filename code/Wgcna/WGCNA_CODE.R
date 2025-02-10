
install.packages("FactoMineR")
install.packages("factoextra")
library(tidyverse)
library(WGCNA)
library(biomaRt)
library(FactoMineR)
library(factoextra)

file_path <- 'D:\\metabolic network\\TPM_GSE126848.txt'
fpkm_data <- read.csv(file_path, sep='\t',row.names = 1)

fpkm_log2 <- fpkm_data %>%
  mutate(across(everything(), ~ log2(. + 1)))


#MAD_data <- fpkm_log2[order(apply(fpkm_log2,1,mad), decreasing = T)[1:5000],]

file_path <- 'D:\\metabolic network\\\Wgcna\\condition.csv'
datTraits <- read.csv(file_path, sep=',',row.names = 1)

datExpr0 <- as.data.frame(t(fpkm_log2))

gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
if(T){

  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1, cex.lab = 1)
  

  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                  colors = rainbow(length(table(datTraits$group))), 
                                  signed = FALSE)

  par(mar = c(1,4,3,1), cex = 0.8)
  

  png("D:\\metabolic network\\\step1_Sample dendrogram and trait_fibrosis_stage_GSE126848.png", 
      width = 16, height = 8, units = "in", res = 300)
  
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait" )
  
  ## Plot a line to show the cut
  # abline(h = 23500, col = "red")
  dev.off()
}



if(T){
  clust <- cutreeStatic(sampleTree, cutHeight = 100, minSize = 10) 
  table(clust)
  keepSamples <- (clust==1)
  datExpr_1 <- datExpr0[keepSamples, ]
  datTraits_1 <- datTraits[keepSamples,]
  dim(datExpr_1) 
}

group_list <- datTraits_1$group
dat.pca <- PCA(datExpr_1, graph = F) 
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point","text"), #"point","text"
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE, 
                    col.ind = group_list, 
                    axes.linetype=NA,  # remove axeslines
                    mean.point=T
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1)
pca
ggsave(filename = "D:\\metabolic network\\\step1_Sample PCA analysis.png", 
       plot = pca, 
       width = 16, 
       height = 8, 
       units = "in", 
       dpi = 300)  
datExpr <-  datExpr_1 
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits_1,file="D:\\metabolic network\\\Wgcna\\step1_input.Rdata")                 





############################### 2.power ###################################
rm(list = ls())  
load("D:\\metabolic network\\\Wgcna\\step1_input.Rdata")
R.sq_cutoff = 0.8 
if(T){
  # Call the network topology analysis function

  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  #SFT.R.sq > 0.8 , slope ≈ -1
  png("D:\\metabolic network\\\step2_power-value.png", 
      width = 16, height = 12, units = "in", res = 300)
  #pdf("D:\\metabolic network\\\step2_power-value.pdf",width = 16,height = 12)
  # Plot the results: 
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # this line corresponds to using an R^2 cut-off of h
  abline(h=R.sq_cutoff ,col="red")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  abline(h=100,col="red")
  dev.off()
}

sft$powerEstimate  #查看估计的最佳power
# power = sft$powerEstimate
power = sft$powerEstimate


if(is.na(power)){

  type = "unsigned"
  nSamples=nrow(datExpr)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}


##################### 3.一步法构建加权共表达网络，识别基因模块 ####################
if(T){
  net <- blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr),
    corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 30,    ##越大模块越少
    mergeCutHeight = 0, ##越大模块越少
    numericLabels = TRUE, 
    saveTOMs = T,
    saveTOMFileBase = "TOM",
    verbose = 3
  )
  table(net$colors) 

}


if(T){
  # Convert labels to colors for plotting
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  # Plot the dendrogram and the module colors underneath
  #pdf("D:\\metabolic network\\\step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
  png("D:\\metabolic network\\\step3_genes-modules_ClusterDendrogram.png", 
      width = 16, height = 12, units = "in", res = 300)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}
save(net, moduleColors, file = "D:\\metabolic network\\\step3_genes_modules.Rdata")


loadedTOM <- load("TOM-block.1.RData")
TOM_from_net <- as.matrix(loadedTOM)



  ### 保存数据
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  # Save module colors and labels for use in subsequent parts
  save(MEs, moduleLabels, moduleColors, geneTree, 
       file = "D:\\metabolic network\\\step3_stepByStep_genes_modules.Rdata")
  
}
load(file = "D:\\metabolic network\\\step3_genes_modules.Rdata")

if(T){
  datTraits_1$group <- as.factor(datTraits_1$group)
  design <- model.matrix(~0+datTraits_1$group)
  colnames(design) <- levels(datTraits_1$group) #get the group
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #Calculate module eigengenes.
  MEs <- orderMEs(MES0)  #Put close eigenvectors next to each other
  moduleTraitCor <- cor(MEs,design,use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf("D:\\metabolic network\\\step4_Module-trait-relationship_heatmap_Stage_NEW1.pdf",
      width = 2*length(colnames(design)), 
      height = 0.6*length(names(MEs)) )
  par(mar=c(5, 9, 3, 3)) #留白：下、左、上、右
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.5,
                 zlim = c(-1,1), 
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = "D:\\metabolic network\\\step4_design.Rdata")
}


if (T) {
  datTraits_1$group <- as.factor(datTraits_1$group)
  design <- model.matrix(~0 + datTraits_1$group)
  colnames(design) <- levels(datTraits_1$group) 
  
  MES0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes 
  MEs <- orderMEs(MES0)  
  

  moduleTraitCor <- cor(MEs, design, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
 
  textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
  dim(textMatrix) <- dim(moduleTraitCor)

  nashIndex <- which(colnames(design) == "NASH")
  moduleTraitCor_NASH <- moduleTraitCor[, nashIndex, drop = FALSE]
  textMatrix_NASH <- textMatrix[, nashIndex, drop = FALSE]
  

  topModules <- order(abs(moduleTraitCor_NASH), decreasing = TRUE)[1:10]
  
 
  moduleTraitCor_top10 <- moduleTraitCor_NASH[topModules, , drop = FALSE]
  textMatrix_top10 <- textMatrix_NASH[topModules, , drop = FALSE]

  svg("D:\\metabolic network\\\step4_Module-trait-relationship_heatmap_NASH_TOP20.svg",
      width = 4,  # 增大图形宽度
      height = 6) # 增大图形高度
  par(mar = c(4, 6, 2, 2))  # 减小边距
  blue_white_orange <- colorRampPalette(c("#7FABD1", "white", "#F27873"))(50)
  labeledHeatmap(Matrix = moduleTraitCor_top10,
                 xLabels = colnames(design)[nashIndex],
                 yLabels = rownames(moduleTraitCor_top10),
                 ySymbols = rownames(moduleTraitCor_top10),
                 colorLabels = F,
                 colors = blue_white_orange,
                 textMatrix = textMatrix_top10,
                 setStdMargins = T,  # 设置标准边距
                 cex.text = 0.5,
                 zlim = c(-1, 1),
                 main = "Top 10 Module-NASH Relationships")
  dev.off()
  
  save(design, file = "D:\\metabolic network\\\step4_design.Rdata")
}




if (T) {

  mes_group <- merge(MEs, datTraits_1, by = "row.names")
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra) 
  
  draw_ggboxplot <- function(data, Module = "Module", group = "group") {
    ggboxplot(data, x = group, y = Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = c("#00AFBB", "#E7B800", "#FC4E07"), # 自定义颜色方案
              legend = "") + stat_compare_means() # 添加统计比较
  }
  

  moduleTraitCor <- cor(MEs, design, use = "p")
  nashIndex <- which(colnames(design) == "NASH")
  moduleTraitCor_NASH <- moduleTraitCor[, nashIndex, drop = FALSE]
  

  topModules <- order(abs(moduleTraitCor_NASH), decreasing = TRUE)[1:10]
  colorNames <- names(MEs)[topModules] # 获取前10个模块的名称
  

  svg("D:\\metabolic network\\\step4_Module-trait-relationship_boxplot_NASH_TOP10.svg", 
      width = 15, 
      height = 1.6 * length(topModules))
  p <- lapply(colorNames, function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange, c(p, ncol = 2)) # 每行排布2个图
  dev.off()
}



if(T){
  mes_group <- merge(MEs,datTraits_1,by="row.names")
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra) 
  draw_ggboxplot <- function(data,Module="Module",group="group"){
    ggboxplot(data,x=group, y=Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = "jco",
              #add="jitter",
              legend = "") +stat_compare_means()
  }

  colorNames <- names(MEs)
  pdf("D:\\metabolic network\\\step4_Module-trait-relationship_boxplot_stage.pdf", width = 15,height = 1.6*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=2)) #排布为每行2个
  dev.off()
}



levels(datTraits_1$group)
choose_group <- "NASH"  

if(T){
  modNames <- substring(names(MEs), 3)
  

  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)


  # trait <- datTraits$groupNo  

  trait <- as.data.frame(design[,choose_group])
  
  geneTraitSignificance <- as.data.frame(cor(datExpr,trait,use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
  names(geneTraitSignificance) <- paste0("GS")
  names(GSPvalue) <- paste0("GS")
  
  #selectModule<-c("blue","green","purple","grey")  
  selectModule <- modNames  
  pdf("D:\\metabolic network\\\step4_gene-Module-trait-significance_stage.pdf",width=7, height=1.5*ncol(MEs))
  par(mfrow=c(ceiling(length(selectModule)/2),2)) #批量作图开始
  for(module in selectModule){
    column <- match(module,selectModule)
    print(module)
    moduleGenes <- moduleColors==module
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for trait",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  }
  dev.off()
}
module = "greenyellow"
gene <- colnames(datExpr) 
inModule <- moduleColors==module
modgene <- gene[inModule]
geneInfo <- data.frame(Gene = gene, Module = moduleColors)
write.csv(geneInfo, file = "D:\\metabolic network\\\Module_Genes_NO_merge_GSE126848.csv", row.names = FALSE)

dissTOM=1-TOM_matrix
save(dissTOM, file = "D:\\metabolic network\\\step5_prepraction.Rdata")
load("D:\\metabolic network\\\step5_prepraction.Rdata")
# 假设 moduleColors 是一个向量，TOM_matrix 是一个矩阵
geneOrder <- net$dendrograms[[1]]$order

# 根据 geneOrder 对 net$colors 进行排序
orderedColors <- net$colors[geneOrder]

# 确保 orderedColors 与 dendrograms 的维度一致
if (length(orderedColors) == length(geneOrder)) {
  print("维度一致，可以继续分析")
} else {
  stop("维度仍然不一致，请检查数据")
}


if(T){
  #TOM=TOMsimilarityFromExpr(datExpr,power=power)
  #dissTOM=1-TOM_matrix
  moduleColors <- labels2colors(orderedColors)
  ## draw all genes 
  if(T){
    geneTree = net$dendrograms[[1]]
    #geneTree = hclust(as.dist(dissTOM), method = "average");
    plotTOM = dissTOM^7
    diag(plotTOM)=NA
    png("D:\\metabolic network\\\step5_TOMplot_Network-heatmap.png",width = 800, height=600) 
    TOMplot(plotTOM,geneTree,moduleColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot")
    dev.off()
  }
  ### draw selected genes to save time...just for test...
  if(F){
    nSelect =0.1*nGenes
    set.seed(123)
    select=sample(nGenes,size = nSelect)
    selectTOM = dissTOM[select,select]
    selectTree = hclust(as.dist(selectTOM),method = "average")
    selectColors = moduleColors[select]
    plotDiss=selectTOM^7
    diag(plotDiss)=NA
    pdf("D:\\metabolic network\\\step5_select_TOMplot_Network-heatmap.pdf",width=8, height=6)
    TOMplot(plotDiss,selectTree,selectColors,
            col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'),
            main="Network heapmap plot of selected gene")
    dev.off()
  }
}


load(net$TOMFiles[1], verbose=T)

## Loading objects:
##   TOM

TOM <- as.matrix(TOM)

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

# 这一部分特别耗时，行列同时做层级聚类
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, all genes")




### 模块相关性展示 Eigengene-adjacency-heatmap
if(T){
  MEs = moduleEigengenes(datExpr_1,moduleColors)$eigengenes
  MET = orderMEs(MEs)
  # 若添加表型数据
  if(T){
    ## 连续型性状
    # MET = orderMEs(cbind(MEs,datTraits$groupNo))
    ## 非连续型性状，需将是否属于这个表型进行0,1数值化，已存于design中
    design
    NASH = as.data.frame(design[,2])
    names(NASH) = "NASH"
    # Add the weight to existing module eigengenes
    MET = orderMEs(cbind(MEs, NASH))
  }
  pdf("D:\\metabolic network\\\step5_module_cor_Eigengene-dendrogram_control.pdf",width = 8,height = 10)
  plotEigengeneNetworks(MET, setLabels="", 
                        marDendro = c(0,4,1,4),  # 留白：下右上左
                        marHeatmap = c(5,5,1,2), # 留白：下右上左
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}



### 条件设置
OrgDb = "org.Mm.eg.db"  # "org.Mm.eg.db"  "org.Hs.eg.db"
genetype = "SYMBOL"    # "SYMBOL"   "ENSEMBL"
table(moduleColors)
choose_module <- c("lightcyan","cyan","greenyellow")

inModule <- moduleColors==module
modgene <- gene[inModule]
