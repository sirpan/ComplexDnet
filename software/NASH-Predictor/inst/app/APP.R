library(shiny)
options(shiny.maxRequestSize = 500*1024^2)  # 500MB
library(WGCNA)
library(data.table)
library(tibble)
library(pheatmap)
library(DT)
library(reticulate)
library(shinyjs)
# 加载自定义函数

source(".//WGCNA//WGCNA_CODE.R")

ui <- fluidPage(
  useShinyjs(),
  theme = shinythemes::shinytheme("cerulean"),
  titlePanel(div(icon("dna"), "NASH-Predictor Analytical Platform", 
                 style = "color: #2C3E50; font-family: 'Arial'; font-size: 28px;")),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      style = "background-color: #F8F9FA; border-right: 2px solid #E9ECEF;",
      
      h4(icon("upload"), "DATA INPUT", 
         style = "color: #3498DB; font-size: 18px; font-weight: bold; margin-top: 15px;"),
      
      # 基因表达数据上传
      fileInput("exprFile", "Gene Expression Data (CSV)",
                accept = ".csv", buttonLabel = "Upload..."),
      
      # 表型数据上传
      fileInput("traitFile", "Phenotypic Data (CSV)",
                accept = ".csv", buttonLabel = "Upload..."),
      
      # 基因集文件上传（带默认值提示）
      fileInput("geneSet", HTML("Gene Set File (CSV)<br/><small>Default: GeneCards_clock_big.csv</small>"), 
                accept = ".csv", buttonLabel = "Upload..."),
      
      # 通路文件上传（带默认值提示）
      fileInput("pathwaySet", HTML("Pathway File (CSV)<br/><small>Default: KEY_PATHWAY.csv</small>"),
                accept = ".csv", buttonLabel = "Upload..."),
      
      hr(style = "border-color: #E9ECEF;"),
      
      h4(icon("sliders"), "ANALYSIS PARAMETERS", 
         style = "color: #3498DB; font-size: 18px; font-weight: bold;"),
      radioButtons("geneNumber", "Gene Selection Strategy:",
                   choices = c("Top 2000 Genes" = 2000,
                               "Top 5000 Genes" = 5000,
                               "Top 10000 Genes" = 10000,
                               "All Genes" = "keep_all"),
                   selected = "keep_all"),
      
      actionButton("run", "Run", 
                   icon = icon("play-circle"), 
                   class = "btn-primary btn-lg",
                   style = "width: 100%; background-color: #18BC9C; border: none; margin-top: 20px;"),
      
      hr(style = "border-color: #E9ECEF;"),
      
      div(icon("info-circle"), "USER GUIDE",
          style = "color: #95A5A6; font-size: 14px; margin-top: 15px; font-weight: bold;"),
      helpText("1. Upload data in CSV format", 
               style = "color: #7B8A8B; font-size: 12px;"),
      helpText("2. Gene selection based on variance ranking",
               style = "color: #7B8A8B; font-size: 12px;"),
      helpText("3. Ensure consistent column headers with example files",
               style = "color: #7B8A8B; font-size: 12px;")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "mainTabs",
        type = "pills",
        
        tabPanel(icon("database"), div("DATA OVERVIEW", 
                                       style = "font-size: 16px; font-weight: bold;"),
                 fluidRow(
                   column(6, 
                          wellPanel(
                            style = "background-color: white;",
                            h4(icon("table"), "Expression Data Preview",
                               style = "font-size: 16px; color: #2C3E50;"),
                            div(DT::DTOutput("dataPreview"), 
                                style = "font-size: 90%; height: 200px; overflow-y: auto;")
                          )),
                   column(6,
                          wellPanel(
                            style = "background-color: white;",
                            h4(icon("clipboard-list"), "Phenotype Data Preview",
                               style = "font-size: 16px; color: #2C3E50;"),
                            div(DT::DTOutput("conditionPreview"),
                                style = "font-size: 90%; height: 200px; overflow-y: auto;")
                          ))
                 ),
                 fluidRow(
                   column(6,
                          wellPanel(
                            style = "background-color: white;",
                            h4(icon("project-diagram"), "Sample Distribution (PCA)",
                               style = "font-size: 16px; color: #2C3E50;"),
                            plotOutput("pcaPlot", height = "400px")
                          )),
                   column(6,
                          wellPanel(
                            style = "background-color: white;",
                            h4(icon("tree"), "Sample Clustering Tree",
                               style = "font-size: 16px; color: #2C3E50;"),
                            plotOutput("sampleCluster", height = "400px")
                          ))
                 )
        ),
        
        tabPanel(icon("sitemap"), div("MODULE ANALYSIS", 
                                      style = "font-size: 16px; font-weight: bold;"),
                 fluidRow(
                   column(6,
                          wellPanel(
                            style = "background-color: white;",
                            h4(icon("fire"), "Module Relationship Heatmap",
                               style = "font-size: 16px; color: #2C3E50;"),
                            plotOutput("moduleHeatmap", height = "400px")
                          )),
                   column(6,
                          wellPanel(
                            style = "background-color: white;",
                            h4(icon("project-diagram"), "Module-Gene Relationships",
                               style = "font-size: 16px; color: #2C3E50;"),
                            plotOutput("moduleDendro", height = "400px")
                          )),
                   column(6,
                          wellPanel(
                            style = "background-color: white;",
                            h4(icon("connectdevelop"), "Module-Trait Association Network",
                               style = "font-size: 16px; color: #2C3E50;"),
                            plotOutput("moduleTraitPlot", height = "400px")
                          )),
                   column(6,
                          wellPanel(
                            style = "background-color: white;",
                            h4(icon("sort-amount-down"), "Module Ranking",
                               style = "font-size: 16px; color: #2C3E50;"),
                            div(DT::DTOutput("moduleTable"),
                                style = "font-size: 90%; height: 400px; overflow-y: auto;"),
                 ))
        )),
        
  
        
        tabPanel(icon("hubspot"), div("PPI NETWORK", 
                                      style = "font-size: 16px; font-weight: bold;"),
                 fluidRow(
                   column(12,
                          wellPanel(
                            style = "background-color: white; padding: 20px;",
                            h4(icon("calculator"), "NASH-Score Results Panel",
                               style = "font-size: 16px; color: #2C3E50; margin-bottom: 15px;"),
                            div(DT::DTOutput("resultsTable"),
                                style = "font-size: 95%; border: 1px solid #ECF0F1; height: 500px;")
                          ))
                 )
        )
      )
    )
  )
)


server <- function(input, output, session) {
  
  data <- reactiveValues(
    pre = NULL,
    wgcna = NULL,
    moduleTrait = NULL,
    ppi = list(
      gapr = NULL,
      finder = NULL,
      mapping = NULL
    )
  )
  
  observeEvent(input$run, {
    req(input$exprFile, input$traitFile)
    
    tryCatch({
      # 清除旧结果
      data$pre <- NULL
      data$wgcna <- NULL
      data$moduleTrait <- NULL
      
      withProgress(message = 'Analysing progress', value = 0, {
        # 数据读取
        incProgress(0.1, detail = "Read data...")
        expr_raw <- fread(input$exprFile$datapath) %>% 
          as.data.frame() %>%
          column_to_rownames(var = names(.)[1])
        
        trait_raw <- fread(input$traitFile$datapath) %>% 
          as.data.frame() %>%
          column_to_rownames(var = names(.)[1])
        
        # 数据校验
        validate(
          need(ncol(expr_raw) >= 10, "A minimum of 10 samples is required"),
          need(nrow(trait_raw) == ncol(expr_raw), "Mismatch between sample size and phenotypic data"),
          need("group" %in% colnames(trait_raw), "data must contain group columns")
        )
        
        # 数据预处理
        incProgress(0.1, detail = "Data preprocessing ...")
        data$pre <- wgcna_preprocess(t(expr_raw), trait_raw,input$geneNumber)
        
        # WGCNA分析
        incProgress(0.2, detail = "network construction ...")
        data$wgcna <- wgcna_analysis(data$pre$datExpr, data$pre$power,data$pre$datTraits)
        
        # 模块-性状关联分析
        incProgress(0.3, detail = "Correlation analysis...")
        design <- model.matrix(~0 + trait_raw$group)

        
        data$moduleTrait <- list(
          cor = data$wgcna$moduleTraitCor_top10,
          pval = data$wgcna$textMatrix_top10
        )
        
        # PPI网络构建
        #incProgress(0.8, detail = "构建PPI网络...")
        #py_run_string("from ppi_module import build_ppi, plot_network")
        #ppi <- py$build_ppi(
          #genes = colnames(data$pre$datExpr),
          #species = "human"
        #)
       # data$ppi <- ppi
      })
      
    
  # 计算NASH分数

      withProgress(message = "Generating a PPI network", value = 0.7, {
        cat("\n--- Starting the PPI Process ---\n")
        
        # 加载Python脚本
        incProgress(0.5, detail = "Initialising the Python environment...")
        use_condaenv("NASH_net", required = TRUE)
        source_python("GAPR/Process_ppi.py")  # 直接调用Python脚本
        
        # 获取目标模块基因
        incProgress(0.7, detail = "Extracting modular genes ...")
        target_modules <- head(data$wgcna$selectModule, 3)
        if(length(target_modules) == 0) stop("No WGCNA module selected")
        
        gene <- colnames(data$pre$datExpr) 
        inModule <- data$wgcna$moduleColors %in% target_modules
        modgene <- gene[inModule]
        
        # 创建临时目录
        temp_dir <- "ppi_network"
        if(!dir.exists(temp_dir)) dir.create(temp_dir)
        
        # 写入基因文件（修正部分）
        gene_path <- file.path(temp_dir, "target_genes.csv")
        gene_df <- data.frame(ensembl_id = modgene)
        fwrite(gene_df, gene_path)
        
        # 设置文件路径（转换为Python兼容路径）
        gene_info_path <- normalizePath(file.path(temp_dir, "Human.GRCh38.p13.annot.tsv"))
        ppi_path <- normalizePath(file.path(temp_dir, "2022PPI_35w.tsv"))
        
        # 验证文件存在
        validate(
          need(file.exists(gene_info_path), "Gene annotation file does not exist"),
          need(file.exists(ppi_path), "PPI network file does not exist"),
          need(file.exists(gene_path), "Target gene file not generated"),
          need(length(modgene) > 0, "The target module has no genes")
        )
        
        # 调用Python函数
        
        py$process_data(
          gene_path = normalizePath(gene_path),
          gene_info_path = gene_info_path,
          ppi_path = ppi_path,
          output_dir = normalizePath(temp_dir)
        )
        
        # 读取生成文件
        
        data$ppi$gapr <- fread(file.path(temp_dir, "GAPR_PPI.csv"))
        data$ppi$finder <- readLines(file.path(temp_dir, "Finder_PPI.txt"))
        data$ppi$mapping <- fread(file.path(temp_dir, "Finder_node_mapping.txt"), 
                                  header = FALSE, sep = " ")
        
        showNotification("PPI network generated successfully!", type = "message")
        if (!is.null(input$geneSet)){
          geneSelect <-read.csv(input$geneSet$datapath)}
        else {geneSelect <-read.csv("./test/key_target_set.csv")}
        
       
        if (!is.null(input$pathwaySet)){
          pathwaySelect <-read.csv(input$pathwaySet$datapath)}
        else {pathwaySelect <-read.csv("./test/KEY_PATHWAY.csv")}
        
        print(geneSelect)
        
        
        # 调用Python脚本
        incProgress(0.8, detail = "Perform GAPR score calculations...")
        gapr_dir <- "GAPR"
        if(!dir.exists(gapr_dir)) dir.create(gapr_dir)
        
        # 写入PPI网络文件
        ppi_path <- file.path(gapr_dir, "GAPR_PPI.txt")
        fwrite(data$ppi$gapr, ppi_path, sep = "\t")
        output_dir="./GAPR/GAPR_Results.csv"
        # 配置Python环境
        source_python(normalizePath("GAPR/GAPR_python.py"))
        
        # 运行GAPR算法
        py$gapr_analysis(
          ppi_path = normalizePath(ppi_path),
        
          output_dir = output_dir
        )
        
        showNotification("GAPR analysis completed!", type = "message")
        incProgress(0.9, detail = "Perform Finder score calculations...")
        # 运行Finder算法
        original_wd <- getwd()
        on.exit({
          setwd(original_wd)
          cat("\n has restored the working directory to:", original_wd, "\n")
        })
        
        # 使用绝对路径操作
        finder_code_dir <- normalizePath("FINDER/code/FINDER_CN")
        setwd(finder_code_dir)
        cat("\n current working directory:", getwd(), "\n")
        # 预编译检查机制
        if (!file.exists("realSolver.c")) {
          incProgress(0.2, detail = "Compiling Cython Extensions...")
          compile_log <- system2("python", "setup.py build_ext -i", stdout = TRUE, stderr = TRUE)
          cat("Compilation Log:\n", paste(compile_log, collapse = "\n"), "\n")
        }
        py_run_string("
import tensorflow as tf

                    ")
        source_python("testReal.py")
        # 预编译检查机制
        if (!file.exists("realSolver.c")) {
          incProgress(0.2, detail = "Compiling Cython Extensions...")
          compile_log <- system2("python", "setup.py build_ext -i", stdout = TRUE, stderr = TRUE)
          cat("Compilation Log:\n", paste(compile_log, collapse = "\n"), "\n")
        }
        
        # 安全写入文件
        ppi_path <- normalizePath("Finder_PPI.txt", mustWork = FALSE)
        writeLines(data$ppi$finder, ppi_path)
        cat("PPI file written:", ppi_path, "\n文件大小:", file.size(ppi_path), "bytes\n")
        py$GetSolution(0.01
        )
        py$EvaluateSolution(0.01, 0)
        
        
        # 安全读取结果
        
        
        showNotification("GAPR analysis completed!", type = "message")
        
        #执行NASH-Score计算
        setwd(original_wd)
        incProgress(0.7, detail = "Perform NASH-Score calculations...")
        source_python("./NASH_Score/NASH_Score.py")
        
        py$cal_nash_score(geneSelect,pathwaySelect)
        data$scores<- read.csv("./NASH_Score/NASH_Score.csv")
       
      })
    }, error = function(e){
      showNotification(paste("Failure to Calculate Score:", e$message), type = "error")
    },finally = {
      # 强制清理TensorFlow资源
      py_run_string("
if 'dqn' in locals():
    dqn.sess.close()
    del dqn
tf.keras.backend.clear_session()
tf.compat.v1.reset_default_graph()
                ")
    }
    )}
  )
  
  
  
  # 数据预览
  # 修改server端的表格渲染函数
  output$dataPreview <- DT::renderDT({
    req(input$exprFile)
    dat <- fread(input$exprFile$datapath)
    DT::datatable(
      head(dat, 10),
      options = list(
        scrollX = TRUE,
        dom = 't',
        pageLength = 5
      ),
      caption = "Preview of the first 5 lines of expression data"
    )
  })
  
  output$conditionPreview <- DT::renderDT({
    req(input$traitFile)
    dat <- fread(input$traitFile$datapath)
    DT::datatable(
      head(dat, 5),
      options = list(
        scrollX = TRUE,
        dom = 't',
        pageLength = 5
      ),
      caption = "Preview of the first 5 rows of phenotype data"
    )
  })
  
  
  # PCA可视化
  output$pcaPlot <- renderPlot({
    req(data$pre)
    data$pre$pcaPlot +
      theme(legend.position = "right") +
      scale_color_brewer(palette = "Set1") +
      labs(title = NULL) +  # 移除标题
      theme(plot.title = element_blank())  # 双重保障移除标题
  })
  
  # 样本聚类树
  output$sampleCluster <- renderPlot({
    req(data$pre)
    plot(data$pre$sampleTree, main = "sample clustering tree", 
         xlab = "", sub = "", cex = 0.7)
    abline(h = 15, col = "red")
  })
  # 模块热图
  # 模块热图
  output$moduleHeatmap <- renderPlot({
    req(data$wgcna)
    heatmap_data <- cor(data$wgcna$MEs, use = "p")
    
    pheatmap(
      heatmap_data,
      color = colorRampPalette(c("navy", "white", "firebrick"))(50),
      clustering_method = "average",
      border_color = NA,
      fontsize_row = 12,  # 增大行字体
      fontsize_col = 12,  # 增大列字体
      main = "Module Eigengene Relationships"
    )
  })
  
  
  # 模块树状图
  
  output$moduleDendro <- renderPlot({
    req(data$wgcna)
    
    # 设置绘图参数
    par(mfrow = c(ceiling(length(data$wgcna$selectModule)/2), 2),
        mar = c(4, 4, 2, 1))
    
    # 循环绘制每个模块
    for(module in data$wgcna$selectModule){
      column <- match(module,data$wgcna$selectModule)
      print(module)
      moduleGenes <- data$wgcna$moduleColors==module
      verboseScatterplot(abs(data$wgcna$geneModuleMembership[moduleGenes, column]),
                         abs(data$wgcna$geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = "Gene significance for trait",
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    }
  }
  )
  
  
  # 模块-性状关联图
  output$moduleTraitPlot <- renderPlot({
    req(data$wgcna)  # 确保数据已加载
    
    validate(
      need(!is.null(data$wgcna$moduleTraitCor_top10), "未找到模块特征数据"),
      need(!is.null(data$wgcna$design), "未找到设计矩阵信息"),
      need(!is.null(data$wgcna$textMatrix_top10), "未找到文本矩阵数据")
    )
    
    labeledHeatmap(
      Matrix = data$wgcna$moduleTraitCor_top10,
      xLabels = colnames(data$wgcna$design)[data$wgcna$nashIndex],
      yLabels = rownames(data$wgcna$moduleTraitCor_top10),
      ySymbols = rownames(data$wgcna$moduleTraitCor_top10),
      colorLabels = FALSE,
      colors = colorRampPalette(c("#7FABD3", "white", "#F27873"))(50),
      textMatrix = data$wgcna$textMatrix_top10,
      setStdMargins = TRUE,
      cex.text = 0.8,  # 增大文本字号
      cex.lab = 1.2,   # 坐标轴标签字号
      cex.main = 1.5,  # 标题字号
      mar = c(8, 10, 3, 5),  # 调整边距 (下,左,上,右)
      zlim = c(-1, 1),
      main = "Top 10 Module-phenotype Relationships"
    )
    
  })
  
  
  # 模块表格
  output$moduleTable <- renderDT({
    req(data$moduleTrait)
    cor_df <- as.data.frame(data$wgcna$moduleTraitCor_top10)
    pval_df <- as.data.frame(data$wgcna$moduleTraitP_top10)
    
    res <- sapply(1:ncol(cor_df), function(i) {
      paste0(round(cor_df[,i], 2), " (", round(pval_df[,i], 4), ")")
    })
    
    datatable(
      as.data.frame(res),
      colnames = colnames(cor_df),
      options = list(
        pageLength = 10, 
        dom = 't',
        scrollX = TRUE
      ),
      caption = "Module-trait correlations (p-values in parentheses)"
    )
  })
   
  output$resultsTable <- DT::renderDT({
  req(data$scores)
  
  # 数据重命名处理
  data_renamed <- data$scores %>%
    select(any_of(c("GeneID", "Symbol", "Finder_score", "sim", "Score_total"))) %>%
    rename(
      `Gene ID` = GeneID,
      `Gene Symbol` = Symbol,
      `Finder Score` = Finder_score,
      `Pathway Match` = sim,
      `Total Score` = Score_total
    )
  
  # 计算目标列的索引（0-based）
  total_score_index <- which(names(data_renamed) == "Total Score") - 1
  
  # 确保索引有效
  validate(need(!is.na(total_score_index), "Error: Total Score column not found"))
  
  # 构建DT参数
  dt <- DT::datatable(
    data_renamed,
    options = list(
      pageLength = 10,
      order = list(list(total_score_index, 'desc')),  # 使用预计算索引
      dom = 'Bfrtip',
      scrollX = TRUE,
      columnDefs = list(
        list(className = 'dt-center', targets = '_all')
      )
    ),
    caption = htmltools::tags$caption(
      style = 'caption-side: top; text-align: center; color: #2C3E50;',
      htmltools::h4("Key targets and composite scores")
    ),
    rownames = FALSE
  ) %>%
    DT::formatRound(
      columns = c("Finder Score", "Pathway Match", "Total Score"),
      digits = 3
    ) %>%
    DT::formatStyle(
      'Total Score',
      backgroundColor = DT::styleInterval(
        c(quantile(data_renamed$`Total Score`, 0.75, na.rm = TRUE)),
        c('#FFFFFF', '#E8F8F5')
      )
    )
  
  dt
})

  
  # PPI网络图
  #output$ppiPlot <- renderPlot({
  #req(data$ppi)
  # py$plot_network(data$ppi)
  #})
  
  # Hub基因表
  #output$hubGenes <- renderDT({
  #req(data$ppi)
  #hub_genes <- data.frame(
  #Gene = names(sort(data$ppi$degree, decreasing = TRUE)[1:10]),
  # Degree = sort(data$ppi$degree, decreasing = TRUE)[1:10]
  #)
  
  #datatable(
  # hub_genes,
  #options = list(pageLength = 5, dom = 't'),
  #caption = "Top 10 Hub基因"
  #)
  
}
shinyApp(ui, server)