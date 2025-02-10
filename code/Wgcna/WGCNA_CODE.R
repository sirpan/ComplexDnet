# 加载必要的包
install.packages("FactoMineR")
install.packages("factoextra")
library(tidyverse)
library(WGCNA)
library(biomaRt)
library(FactoMineR)
library(factoextra)
# 读取FPKM数据
file_path <- 'D:\\metabolic network\\转录组\\Wgcna\\TPM_GSE126848.txt'
fpkm_data <- read.csv(file_path, sep='\t',row.names = 1)
# log2转换数据
fpkm_log2 <- fpkm_data %>%
  mutate(across(everything(), ~ log2(. + 1)))


#MAD_data <- fpkm_log2[order(apply(fpkm_log2,1,mad), decreasing = T)[1:5000],]
#读取分组情况
file_path <- 'D:\\metabolic network\\转录组\\Wgcna\\condition.csv'
datTraits <- read.csv(file_path, sep=',',row.names = 1)

datExpr0 <- as.data.frame(t(fpkm_log2))
### 判断数据质量--缺失值
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
  # 针对样本做聚类树
  sampleTree <- hclust(dist(datExpr0), method = "average")
  par(mar = c(0,5,2,0))
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, 
       cex.axis = 1, cex.main = 1, cex.lab = 1)
  
  ## 若样本有性状、表型，可以添加对应颜色，查看是否聚类合理
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                  colors = rainbow(length(table(datTraits$group))), 
                                  signed = FALSE)
  
  ## 绘制样品的系统聚类树及对应性状
  par(mar = c(1,4,3,1), cex = 0.8)
  
  # 使用适当的宽度和高度设置图形设备
  png("D:\\metabolic network\\转录组\\step1_Sample dendrogram and trait_fibrosis_stage_GSE126848.png", 
      width = 16, height = 8, units = "in", res = 300)
  
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      main = "Sample dendrogram and trait" )
  
  ## Plot a line to show the cut
  # abline(h = 23500, col = "red") # 根据实际情况而定
  dev.off()
}


##若存在显著离群点；剔除掉
if(T){
  clust <- cutreeStatic(sampleTree, cutHeight = 100, minSize = 10) # cutHeight根据实际情况而定
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
                    repel = TRUE, #标签不重叠
                    col.ind = group_list, # 分组上色
                    axes.linetype=NA,  # remove axeslines
                    mean.point=T#去除分组中心点
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1) #坐标轴的纵横比
pca
ggsave(filename = "D:\\metabolic network\\转录组\\step1_Sample PCA analysis.png", 
       plot = pca,  # 使用pca对象
       width = 16, 
       height = 8, 
       units = "in",  # 以英寸为单位
       dpi = 300)  # 设置分辨率为300DPI
datExpr <-  datExpr_1 
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
save(nGenes,nSamples,datExpr,datTraits_1,file="D:\\metabolic network\\转录组\\Wgcna\\step1_input.Rdata")                 





############################### 2.挑选最佳阈值power ###################################
rm(list = ls())  
load("D:\\metabolic network\\转录组\\Wgcna\\step1_input.Rdata")
R.sq_cutoff = 0.8  #设置R^2 cut-off，默认为0.85
if(T){
  # Call the network topology analysis function
  #设置power参数选择范围
  powers <- c(seq(1,20,by = 1), seq(22,30,by = 2)) 
  sft <- pickSoftThreshold(datExpr, 
                           networkType = "unsigned",
                           powerVector = powers, 
                           RsquaredCut = R.sq_cutoff,  
                           verbose = 5)
  #SFT.R.sq > 0.8 , slope ≈ -1
  png("D:\\metabolic network\\转录组\\step2_power-value.png", 
      width = 16, height = 12, units = "in", res = 300)
  #pdf("D:\\metabolic network\\转录组\\step2_power-value.pdf",width = 16,height = 12)
  # Plot the results: 寻找拐点，确认最终power取值
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

# 若无向网络在power小于15或有向网络power小于30内，没有一个power值使
# 无标度网络图谱结构R^2达到0.8且平均连接度在100以下，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if(is.na(power)){
  # 官方推荐 "signed" 或 "signed hybrid"
  # 为与原文档一致，故未修改
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
  # power: 上一步计算的软阈值
  # maxBlockSize:计算机能处理的最大模块的基因数量(默认5000),16G内存可以处理2万个，
  # 计算资源允许的情况下最好放在一个block里面。
  # corType：计算相关性的方法；可选pearson(默认)，bicor。后者更能考虑离群点的影响。
  # networkType:计算邻接矩阵时，是否考虑正负相关性；默认为"unsigned",可选"signed", "signed hybrid"
  # TOMType：计算TOM矩阵时，是否考虑正负相关性；默认为"signed",可选"unsigned"。但是根据幂律转换的邻接矩阵(权重)的非负性，所以认为这里选择"signed"也没有太多的意义。
  # numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
  # saveTOMs：最耗费时间的计算，可存储起来供后续使用，
  # mergeCutHeight: 合并模块的阈值，越大模块越少,一般为0.25
  # minModuleSize: 每个模块里最少放多少个基因，设定越大模块越少
  # 输出结果根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
  # **0 (grey)**表示**未**分入任何模块的基因。
}

## 模块可视化，层级聚类树展示各个模块
if(T){
  # Convert labels to colors for plotting
  moduleColors <- labels2colors(net$colors)
  table(moduleColors)
  # Plot the dendrogram and the module colors underneath
  #pdf("D:\\metabolic network\\转录组\\step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
  png("D:\\metabolic network\\转录组\\step3_genes-modules_ClusterDendrogram.png", 
      width = 16, height = 12, units = "in", res = 300)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
}
save(net, moduleColors, file = "D:\\metabolic network\\转录组\\step3_genes_modules.Rdata")


loadedTOM <- load("TOM-block.1.RData")
TOM_from_net <- as.matrix(loadedTOM)

#####################  分布法完成网络构建，一般不用 #################################
if(T){
  ## 构建加权共表达网络分为两步：
  ## 1. 计算邻近值，也是就是两个基因在不同样品中表达量的表达相关系数(pearson correlation rho)，
  ## 2. 计算topology overlap similarity (TOM)。用TOM表示两个基因在网络结构上的相似性，即两个基因如果具有相似的邻近基因，这两个基因更倾向于有相互作用。
  
  ###(1)网络构建 Co-expression similarity and adjacency 
  adjacency = adjacency(datExpr_1, power = power) 
  
  ###(2) 邻近矩阵到拓扑矩阵的转换，Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1-TOM
  
  ###(3) 聚类拓扑矩阵 Clustering using TOM
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  # Plot the resulting clustering tree (dendrogram)
  ## 这个时候的geneTree与一步法的 net$dendrograms[[1]] 性质类似，但是还需要进行进一步处理
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04)
  
  ###(4) 聚类分支的修整 dynamicTreeCut 
  ################# set the minimum module size ##############################
  minModuleSize = 100
  ####
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize)
  table(dynamicMods)
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  
  ###(5) 聚类结果相似模块的融合 Merging of modules whose expression profiles are very similar
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  #一般选择 height cut 为0.25,对应于有75%相关性，进行融合
  ###################### set  Merging height cut  ################################
  MEDissThres = 0.5
  ####
  # Plot the result
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  # 统计mergedmodule
  table(mergedColors)
  
  ### (6) plot the gene dendrogram 
  pdf(file = "D:\\metabolic network\\转录组\\step3_stepbystep_DynamicTreeCut_genes-modules.pdf", width = 16,height = 12)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  ### 保存数据
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder)-1
  MEs = mergedMEs
  # Save module colors and labels for use in subsequent parts
  save(MEs, moduleLabels, moduleColors, geneTree, 
       file = "D:\\metabolic network\\转录组\\step3_stepByStep_genes_modules.Rdata")
  
}
load(file = "D:\\metabolic network\\转录组\\step3_genes_modules.Rdata")
## 模块与表型的相关性热图
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
  
  pdf("D:\\metabolic network\\转录组\\step4_Module-trait-relationship_heatmap_Stage_NEW1.pdf",
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
  save(design, file = "D:\\metabolic network\\转录组\\step4_design.Rdata")
}


if (T) {
  datTraits_1$group <- as.factor(datTraits_1$group)
  design <- model.matrix(~0 + datTraits_1$group)
  colnames(design) <- levels(datTraits_1$group)  # 获取组别信息
  
  MES0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes  # 计算模块特征基因
  MEs <- orderMEs(MES0)  # 将相近的特征基因排在一起
  
  # 计算模块与表型的相关性和P值
  moduleTraitCor <- cor(MEs, design, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # 将相关性和P值组合成文本矩阵
  textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  # 提取与 NASH 相关的列
  nashIndex <- which(colnames(design) == "NASH")
  moduleTraitCor_NASH <- moduleTraitCor[, nashIndex, drop = FALSE]
  textMatrix_NASH <- textMatrix[, nashIndex, drop = FALSE]
  
  # 对模块按与 NASH 的相关性绝对值排序，选择排名前10的模块
  topModules <- order(abs(moduleTraitCor_NASH), decreasing = TRUE)[1:10]
  
  # 选择排名前10的模块的相关性数据
  moduleTraitCor_top10 <- moduleTraitCor_NASH[topModules, , drop = FALSE]
  textMatrix_top10 <- textMatrix_NASH[topModules, , drop = FALSE]
  
  # 调整图形设备的尺寸
  svg("D:\\metabolic network\\转录组\\step4_Module-trait-relationship_heatmap_NASH_TOP20.svg",
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
  
  save(design, file = "D:\\metabolic network\\转录组\\step4_design.Rdata")
}




if (T) {
  # 将模块特征基因与表型数据合并
  mes_group <- merge(MEs, datTraits_1, by = "row.names")
  library(gplots)
  library(ggpubr)
  library(grid)
  library(gridExtra) 
  
  # 定义绘制箱线图的函数，修改颜色方案
  draw_ggboxplot <- function(data, Module = "Module", group = "group") {
    ggboxplot(data, x = group, y = Module,
              ylab = paste0(Module),
              xlab = group,
              fill = group,
              palette = c("#00AFBB", "#E7B800", "#FC4E07"), # 自定义颜色方案
              legend = "") + stat_compare_means() # 添加统计比较
  }
  
  # 计算模块与表型的相关性，用于选择前10个模块
  moduleTraitCor <- cor(MEs, design, use = "p")
  nashIndex <- which(colnames(design) == "NASH")
  moduleTraitCor_NASH <- moduleTraitCor[, nashIndex, drop = FALSE]
  
  # 排序并选择与 NASH 相关性绝对值最高的前10个模块
  topModules <- order(abs(moduleTraitCor_NASH), decreasing = TRUE)[1:10]
  colorNames <- names(MEs)[topModules] # 获取前10个模块的名称
  
  # 输出前10模块的箱线图
  svg("D:\\metabolic network\\转录组\\step4_Module-trait-relationship_boxplot_NASH_TOP10.svg", 
      width = 15, 
      height = 1.6 * length(topModules))
  p <- lapply(colorNames, function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange, c(p, ncol = 2)) # 每行排布2个图
  dev.off()
}


### 模块与表型的相关性boxplot图 
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
  # 批量画boxplot
  colorNames <- names(MEs)
  pdf("D:\\metabolic network\\转录组\\step4_Module-trait-relationship_boxplot_stage.pdf", width = 15,height = 1.6*ncol(MEs))
  p <- lapply(colorNames,function(x) {
    draw_ggboxplot(mes_group, Module = x, group = "group")
  })
  do.call(grid.arrange,c(p,ncol=2)) #排布为每行2个
  dev.off()
}


### 基因与模块、表型的相关性散点图
#所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因算出相关系数， 
#如果跟性状显著相关的基因也跟某个模块显著相关，那么这些基因可能就非常重要。

# 选择离散性状的表型
levels(datTraits_1$group)
choose_group <- "NASH"  

if(T){
  modNames <- substring(names(MEs), 3)
  
  ### 计算模块与基因的相关性矩阵 
  ## Module Membership: 模块内基因表达与模块特征值的相关性
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste0("MM", modNames)
  names(MMPvalue) <- paste0("p.MM", modNames)
  
  ###  计算性状与基因的相关性矩阵 
  ## Gene significance，GS：比较样本某个基因与对应表型的相关性
  ## 连续型性状
  # trait <- datTraits$groupNo  
  ## 非连续型性状，需转为0-1矩阵, 已存于design中
  trait <- as.data.frame(design[,choose_group])
  
  geneTraitSignificance <- as.data.frame(cor(datExpr,trait,use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance),nSamples))
  names(geneTraitSignificance) <- paste0("GS")
  names(GSPvalue) <- paste0("GS")
  
  ### 可视化基因与模块、表型的相关性.
  #selectModule<-c("blue","green","purple","grey")  ##可以选择自己想要的模块
  selectModule <- modNames  ## 全部模块批量作图
  pdf("D:\\metabolic network\\转录组\\step4_gene-Module-trait-significance_stage.pdf",width=7, height=1.5*ncol(MEs))
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
write.csv(geneInfo, file = "D:\\metabolic network\\转录组\\Module_Genes_NO_merge_GSE126848.csv", row.names = FALSE)

dissTOM=1-TOM_matrix
save(dissTOM, file = "D:\\metabolic network\\转录组\\step5_prepraction.Rdata")
load("D:\\metabolic network\\转录组\\step5_prepraction.Rdata")
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
    png("D:\\metabolic network\\转录组\\step5_TOMplot_Network-heatmap.png",width = 800, height=600) 
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
    pdf("D:\\metabolic network\\转录组\\step5_select_TOMplot_Network-heatmap.pdf",width=8, height=6)
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
  pdf("D:\\metabolic network\\转录组\\step5_module_cor_Eigengene-dendrogram_control.pdf",width = 8,height = 10)
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
