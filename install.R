#' 环境健康检查与智能安装系统
initializeComplexDnet <- function() {
  # 系统级依赖检测
  if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman")
  }
  
  # 配置并行加载与错误恢复机制
  suppressWarnings({
    # CRAN源配置
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    
    # 检测并安装BiocManager
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
      if (!requireNamespace("BiocManager")) stop("BiocManager安装失败，请检查网络连接！")
    }
    
    # 安装关键R包（带版本检测）
    pkgStatus <- pacman::p_loaded()
    requiredCRAN <- c("devtools", "shiny", "shinythemes")
    requiredBioc <- c("GO.db", "preprocessCore", "impute", "biomaRt")
    
    # 智能批量安装（不重复安装）
    if (!all(c(requiredCRAN, requiredBioc) %in% pkgStatus)) {
      pacman::p_load(
        char = requiredCRAN,
        install = TRUE,
        update = FALSE
      )
      BiocManager::install(requiredBioc, ask = FALSE, update = FALSE)
    }
    
    # Python环境构建器
    configurePythonEnv <- function() {
      if (!requireNamespace("reticulate")) install.packages("reticulate")
      
      envName <- "ComplexDnet-env"
      if (!envName %in% reticulate::virtualenv_list()) {
        message("正在构建隔离Python环境...")
        reticulate::virtualenv_create(envName, python_version = "3.8")
      }
      
      reticulate::use_virtualenv(envName, required = TRUE)
      # 安装PyPi依赖
      reqPyPkg <- c("scanpy==1.9.1", "anndata==0.8.0")
      installed <- reticulate::py_list_packages(envName)
      if (!all(reqPyPkg %in% installed$package)) {
        reticulate::py_install(reqPyPkg, envname = envName) 
      }
    }
    
    # 条件执行Python配置
    if (interactive()) tryCatch(configurePythonEnv(),
                                error = function(e) message("Python配置跳过：", e))
  })
  
  # 动态加载主功能模块
  if (!require("ComplexDnet")) {
    message("正在从GitHub安装ComplexDnet...")
    devtools::document()
  }
}

#' 应用启动器增强版
launchApp <- function() {
  if (!require("ComplexDnet")) {
    stop("主程序包加载失败，请检查安装日志！")
  }
  
  # 环境验证检查
  dependencyCheck <- ComplexDnet::validateDependencies()
  if (!all(dependencyCheck$status)) {
    warning("关键依赖缺失：", 
            paste(names(dependencyCheck$status[!dependencyCheck$status]), collapse = ", "))
  }
  
  # 带错误恢复的应用启动
  tryCatch(
    ComplexDnet::runComplexDnetApp(launch.browser = TRUE),
    error = function(e) {
      message("应用启动异常：", e)
      message("尝试清洁安装：remotes::install_github('yourusername/ComplexDnet', force=TRUE)")
    }
  )
}

# 执行主流程
if (interactive()) {
  initializeComplexDnet()
  launchApp()
}
