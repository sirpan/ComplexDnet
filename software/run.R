#Install some dependency packages

install.packages("BiocManager")
BiocManager::install(c("GO.db", "preprocessCore", "impute"))

BiocManager::install("biomaRt")

install.packages("shinythemes")
BiocManager::install(c("GO.db", "preprocessCore", "impute"))

BiocManager::install("biomaRt")

install.packages("shinythemes")

# Installation of the necessary tools
install.packages(c("devtools", "roxygen2"))
# Configuring or Removing Associated Python Environments
# if ("NASH_net" %in% reticulate::conda_list()$name) {
#   reticulate::conda_remove("NASH_net")
# }
# Set the working directory to the current NASHpredictor root directory
setwd("D:\\APP\\NASHpredictor")

# Generating documents and NAMESPACE
devtools::document()

# Formal packaging (automatic processing of paths)
devtools::build(path = "D:/APP", binary = TRUE) 

# Load Test
library(NASHpredictor)
  
runNASHApp()         
