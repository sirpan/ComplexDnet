library(shiny)
library(WGCNA)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(biomaRt)
library(tidyverse)
#--------------------- 功能函数封装 ---------------------#
wgcna_preprocess <- function(datExpr, datTraits,num) {
  progress <- Progress$new()

  # log2转换数据
  
  datExpr <- log2(datExpr + 1)
  # 基因过滤
  gsg <- goodSamplesGenes(datExpr,verbose = 3)
  gsg$allOK
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes],
                                                collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:",
                       paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  #添加基因变异度过滤（保留高变异基因）
  gene_vars <- apply(datExpr, 2, var)
  if (num == "keepall") {
    datExpr <- datExpr
  } else {
    keep <- order(gene_vars, decreasing = TRUE)[1:num]
    datExpr <- datExpr[, keep]
  }
  
  
  # 样本聚类
  sampleTree <- hclust(dist(datExpr), method = "average")
  clust <- cutreeStatic(sampleTree, cutHeight = 100, minSize = 10) # cutHeight根据实际情况而定
  keepSamples <- (clust==1)
  datExpr <- datExpr[keepSamples, ]
  datTraits <- datTraits[keepSamples,]


  
  # PCA分析
  pca_plot <- fviz_pca_ind(
    PCA(datExpr, graph = FALSE),
    geom.ind = c("point"),
    col.ind = as.factor(datTraits$group),
    title = "样本PCA分析"
  )
  #挑选阈值
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
  }
  
  sft$powerEstimate  #查看估计的最佳power
  # power = sft$powerEstimate
  power = sft$powerEstimate
  progress$close()
  
  return(list(
    datExpr = datExpr,
    sampleTree = sampleTree,
    pcaPlot = pca_plot,
    power=power,
    datTraits=datTraits
  ))
}

wgcna_analysis <- function(datExpr, power, datTraits, minModuleSize = 30) {

  # 动态参数配置
  wgcna_version <- packageVersion("WGCNA")
  if (wgcna_version >= "1.72") {
    corArgs <- list(
      corFn = WGCNA::cor, 
      corOptions = list(use = 'p', method = 'pearson')
    )
  } else {
    corArgs <- list(corType = 'pearson')
  }
  
  # 网络构建
  net <- do.call(blockwiseModules, c(
    list(
      datExpr = datExpr,
      power = power,
      networkType = "unsigned",
      TOMType = "unsigned",
      minModuleSize = minModuleSize,
      maxBlockSize = ncol(datExpr),
      numericLabels = TRUE,
      saveTOMs = TRUE,
      saveTOMFileBase = "TOM",
      verbose = 3
    ),
    corArgs  # 动态参数
  ))
  
  

  

  
  # 模块特征分析
  MEs <- moduleEigengenes(datExpr, net$colors)$eigengenes
  moduleColors <- labels2colors(net$colors)
  #排名靠前模块提取
  datTraits$group <- as.factor(datTraits$group)
  design <- model.matrix(~0 + datTraits$group)
  colnames(design) <- levels(datTraits$group)  # 获取组别信息
  
  MES0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes  # 计算模块特征基因
  MEs <- orderMEs(MES0)  # 将相近的特征基因排在一起
  nSamples=nrow(datExpr)
  # 计算模块与表型的相关性和P值
  moduleTraitCor <- cor(MEs, design, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # 将相关性和P值组合成文本矩阵
  textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  # 提取与 NASH 相关的列
  nashIndex <- which(colnames(design) == "Disease")
  moduleTraitCor_NASH <- moduleTraitCor[, nashIndex, drop = FALSE]
  moduleTraitPvalue<- moduleTraitPvalue[, nashIndex, drop = FALSE]
  textMatrix_NASH <- textMatrix[, nashIndex, drop = FALSE]
  
  # 对模块按与 NASH 的相关性绝对值排序，选择排名前10的模块
  # 安全选择模块数量
  n_mods <- length(moduleTraitCor_NASH)
  n_select <- min(10, n_mods)
  
  # 按绝对值排序并选择模块
  topModules <- order(abs(moduleTraitCor_NASH), decreasing = TRUE)[seq_len(n_select)]
  # 选择排名前10的模块的相关性数据
  moduleTraitCor_top10 <- moduleTraitCor_NASH[topModules, , drop = FALSE]
  moduleTraitP_top10<- moduleTraitPvalue[topModules, , drop = FALSE]
  textMatrix_top10 <- textMatrix_NASH[topModules, , drop = FALSE]
  #基因和模块的关联
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  
  choose_group <-"Disease" 
  trait <- as.data.frame(design[,choose_group])
  geneTraitSignificance <- as.data.frame(cor(datExpr, trait, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  rownames(moduleTraitCor) <- gsub("^ME", "", rownames(moduleTraitCor))
  # 获取前10个模块
  module_order <- order(-abs(moduleTraitCor[,1]))
  selectModule <- rownames(moduleTraitCor)[module_order[1:4]]

  
  
  return(list(
    net = net,
    MEs = MEs,
    moduleColors = moduleColors,
    moduleTraitCor_top10=moduleTraitCor_top10,
    textMatrix_top10=textMatrix_top10,
    moduleTraitP_top10=moduleTraitP_top10,
    moduleTraitCor=moduleTraitCor,
    moduleTraitPvalue=moduleTraitPvalue,
    design=design,
    nashIndex=nashIndex,
    selectModule=selectModule,
    geneTraitSignificance=geneTraitSignificance,
    geneModuleMembership=geneModuleMembership

    
  ))
}

#--------------------- Shiny应用逻辑 ---------------------#
ui <- fluidPage(
  titlePanel("NASH预测分析系统"),
  sidebarLayout(
    sidebarPanel(
      fileInput("exprFile", "上传表达矩阵(CSV)", accept = ".csv"),
      fileInput("traitFile", "上传表型数据(CSV)", accept = ".csv"),
      numericInput("power", "软阈值", value = 6, min = 1, max = 20),
      actionButton("analyze", "开始分析")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("数据质量", 
                 plotOutput("pcaPlot"),
                 plotOutput("sampleDendro")),
        tabPanel("网络分析",
                 plotOutput("moduleDendro"),
                 plotOutput("moduleHeatmap")),
        tabPanel("关键模块",
                 plotOutput("moduleTraitPlot"),
                 tableOutput("keyModules"))
      )
    )
  )
)

server <- function(input, output) {
  data <- reactiveValues()
  
  observeEvent(input$analyze, {
    req(input$exprFile, input$traitFile)
    
    tryCatch({
      # 数据读取
      expr_data <- read.csv(input$exprFile$datapath, row.names = 1)
      trait_data <- read.csv(input$traitFile$datapath, row.names = 1)
      
      # 数据预处理
      pre <- wgcna_preprocess(t(expr_data), trait_data)
      data$pca <- pre$pcaPlot
      data$sampleTree <- pre$sampleTree
      
      # WGCNA分析
      wgcna_res <- wgcna_analysis(pre$datExpr, input$power)
      data$modules <- wgcna_res$moduleColors
      data$MEs <- wgcna_res$MEs
      
      # 模块-性状关联分析
      design <- model.matrix(~0 + trait_data$group)
      moduleTraitCor <- cor(wgcna_res$MEs, design, use = "p")
      data$keyModules <- head(sort(moduleTraitCor[,1], decreasing = TRUE), 10)
      
    }, error = function(e) {
      showNotification(paste("错误:", e$message), type = "error")
    })
  })
  
  # 可视化输出
  output$pcaPlot <- renderPlot(data$pca)
  
  output$sampleDendro <- renderPlot({
    plot(data$sampleTree, main = "样本聚类树", xlab = "")
  })
  
  output$moduleDendro <- renderPlot({
    plotDendroAndColors(
      data$modules$dendrograms[[1]],
      data$modules$colors,
      "Module Colors",
      dendroLabels = FALSE
    )
  })
  
  output$moduleHeatmap <- renderPlot({
    labeledHeatmap(
      Matrix = cor(data$MEs, use = "p"),
      xLabels = names(data$MEs),
      yLabels = names(data$MEs),
      colorLabels = FALSE,
      colors = blueWhiteRed(50),
      main = "模块相关性热图"
    )
  })
  
  output$keyModules <- renderTable({
    data.frame(
      模块 = names(data$keyModules),
      相关性 = round(data$keyModules, 3)
    )
  }, rownames = FALSE)
}

shinyApp(ui, server)
