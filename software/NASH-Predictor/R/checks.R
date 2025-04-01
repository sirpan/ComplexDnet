#' @export
# R/runApp.R
runNASHApp <- function() {
  # 环境检测增强版
  if (!"NASH_net" %in% reticulate::conda_list()$name) {
    if (interactive()) {
      ans <- utils::askYesNo("检测到未安装Python环境，现在安装依赖？")
    } else {
      ans <- FALSE  # 非交互环境自动跳过
    }
    
    if (isTRUE(ans)) {
      install_nash_python_deps()
    } else {
      message("操作取消，请手动执行 install_nash_python_deps() 安装环境")
      return(invisible(FALSE))
    }
  }
  
  # 环境激活验证
  tryCatch({
    reticulate::use_condaenv("NASH_net", required = TRUE)
    message("当前Python环境：", reticulate::py_config()$python)
  }, error = function(e) {
    stop("环境激活失败：", e$message)
  })
  
  # 启动Shiny应用
  app_dir <- system.file("app", package = "NASHpredictor")
  if (app_dir == "") {
    stop("无法找到Shiny应用目录")
  }
  shiny::runApp(app_dir, launch.browser = TRUE)
}
