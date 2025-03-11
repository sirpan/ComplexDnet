# 采用CIBERSORT中处理的LUAD的FPKM数据作为基因表达矩阵
load("D:\\metabolic network\\转录组\\Wgcna\\step1_input.Rdata")
##        TCGA-35-5375-01A-01R-1628-07 TCGA-55-A4DF-01A-11R-A24H-07 TCGA-95-8039-01A-11R-2241-07 TCGA-MP-A4T4-01A-11R-A262-07 TCGA-62-A471-01A-12R-A24H-07
## A2M                         40.9806                      45.5385                     226.9787                     125.4271                      23.3876
## A4GALT                       6.3743                       7.8681                       2.7125                       9.7885                      20.6624
## AAAS                         6.3150                       8.2626                       8.2346                       8.3154                      12.2582
## AACS                         2.6022                       3.5811                       1.5165                       2.0158                       5.4609
## AADAC                       28.5316                       1.8770                      10.4736                       0.0450                       4.0560

# 从TISIDB下载免疫基因集列表（http://cis.hku.hk/TISIDB/data/download/CellReports.txt）
library(tidyverse)

cellMarker <- read.csv("D:\\metabolic network\\转录组\\Wgcna\\CellReports.txt", header = F, sep = "\t") # 用EXCEL打开删除NA列
cellMarker <- cellMarker %>% column_to_rownames("V1") %>% t()

a <- cellMarker
a <- a[1:nrow(a), ]
set <- colnames(a)
geneSet <- list()
# i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  geneSet[[i]] <-x
}
BiocManager::install("GSVA")
library(GSVA)
BiocManager::install("limma")
BiocManager::install("GSEABase")
library(limma)
library(GSEABase)
datExpr_t=t(datExpr)
ssgsea <- gsva(datExpr_t, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
# 创建设计矩阵
design <- model.matrix(~ datTraits_1$group)

# 线性模型分析
fit <- lmFit(ssgsea, design)
fit <- eBayes(fit)

# 提取显著性结果
results <- topTable(fit, coef = 2, number = nrow(ssgsea), adjust.method = "fdr")
write.csv(results, file = "D:\\metabolic network\\转录组\\wgcna\\ssGSEA_P_VALUE.csv", row.names = TRUE)
# 筛选出P值小于0.05的细胞类型
significant_results <- results[results$adj.P.Val < 0.05, ]
a <- ssgsea %>% t() %>% as.data.frame()
a$group <- datTraits_1$group
a <- a %>% rownames_to_column("sample")

write.table(a, "D:\\metabolic network\\转录组\\Wgcna\\ssGSEA.txt", sep = "\t", row.names = T, col.names = NA, quote = F)

# Min-Max标准化是指对原始数据进行线性变换，将值映射到[0，1]之间
# 这里是将每个样本中不同的免疫细胞比例标准化到0-1之间
ssgsea.1 <- ssgsea
for (i in colnames(ssgsea)) {
  #i <- colnames(ssgsea)[1]
  ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))
  
}
apply(ssgsea.1[,1:6], 2, range)
library(ggpubr)


LUAD_ggsea <- gather(a,key = ssgsea, value = Expression, -c(group,sample)) 

ggplot(LUAD_ggsea, aes(x = ssgsea, y = Expression)) + 
  labs(y="Expression", x =  NULL) +  
  geom_boxplot(aes(fill = group), position = position_dodge(0.5), width = 0.5, outlier.alpha = 0) + 
  scale_fill_manual(values = c("#096EA9", "#B33D27")) +
  theme_bw() + 
  theme(plot.title = element_text(size = 12,color="black",hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1 ),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.text = element_text(size= 12),
        legend.title= element_text(size= 12)) + 
  stat_compare_means(aes(group =  group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)

# 基因和免疫细胞相关性分析
hub_gene <- read.csv("D:\\metabolic network\\转录组\\Wgcna\\hub_gene_exp.csv",row.names = 1)
hub_gene=t(hub_gene)
ssgsea_t=t(ssgsea.1)
correlation_matrix <- cor(hub_gene, ssgsea_t, method = "spearman")
write.table(correlation_matrix, "D:\\metabolic network\\转录组\\Wgcna\\ssGSEA_hub_gene_correlation.txt", sep = "\t", row.names = T, col.names = NA, quote = F)
