# R/install_python.R
#' 安装Python依赖环境
#' @export
install_nash_python_deps <- function() {

  envs <- reticulate::conda_list()
  if (!"NASH_net" %in% envs$name) {
    message("Creating a Python 3.6 environment...")
    reticulate::conda_create(envname = "NASH_net", python_version = "3.6")
  }
  # 安装依赖包
  message("Configuring Tsinghua Mirror Source Accelerated Installation...")
  req_file <- system.file("python", "requirements.txt", package = "NASHpredictor")
  if (!file.exists(req_file)) {
    stop("The dependency file does not exist:", req_file, "\n",
         "Make sure the package is installed correctly and includes inst/python/requirements.txt.")
  }
  
  message("Start installing Python dependency packages...(may take 5-10 minutes)")
  tryCatch(
    reticulate::py_install(
      envname = "NASH_net",
      packages = c("-r", req_file),
      pip = TRUE,
      pip_options = c(
        "--default-timeout=1000",
        "-i", "https://pypi.tuna.tsinghua.edu.cn/simple",
        "--trusted-host pypi.tuna.tsinghua.edu.cn"
      )
    ),
    error = function(e) {
      stop("Dependency installation failed:", e$message, "\n",
           "Please check your network connection and retry, or install manually:pip install -r ", req_file)
    }
  )

  # 验证安装
  message("\n验证环境配置...")
  tryCatch({
    reticulate::use_condaenv(env_name, required = TRUE)
    check_python_deps()  # 假设存在验证函数
    message("\033[32mSuccess! 环境验证通过\033[39m")
    return(invisible(TRUE))
  }, error = function(e) {
    message("\033[31m验证失败: ", e$message, "\033[39m")
    return(invisible(FALSE))
  })
}