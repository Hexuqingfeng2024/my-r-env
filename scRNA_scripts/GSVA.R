###GSVA分析
# H：包含了由多个已知的基因集构成的超基因集，每个H 类别的基因集都对应多个基础的其他类别的基因集
# C1：包含人类每条染色体上的不同cytoband区域对应的基因集合。根据不染色体编码进行二级分类
# C2：包含了在线通路数据库、PubMed 出版物和领域专家知识的精选基因集
# C3：基于对 miRNA 种子序列和预测的转录因子结合位点的基因靶点预测的调控靶基因组
# C4：包含了计算机软件预测出来的基因集合，主要是和癌症相关的基因
# C5：包含了由相同GO术语注释的基因集
# C6：包含了致癌特征基因集：直接从来自癌症基因扰动的微阵列基因表达数据中定义
# C7：包含了免疫特征基因集：代表免疫系统内的细胞状态和扰动
# C8：包含了细胞类型特征基因组：从人体组织单细胞测序研究中确定的簇标记中收集

library(GSVA)
#BiocManager::install('GSVA')
library(GSEABase)
library(msigdbr)
#①基于单细胞表达量矩阵的gsea和gsva（需要生物学功能数据库）
#②基于细胞亚型表达量矩阵的gsea和gsva（需要生物学功能数据库）

####这里是第①种####
## 取出标准化后的矩阵
genesets <- msigdbr(species = "Homo sapiens",category ='H') 
#genesets <- msigdbr(species = "Homo sapiens") 
table(sce.all$celltype)
sce = subset(sce.all,downsample=200)
table(genesets$gs_cat)
genesets = genesets[genesets$gs_cat %in% c('H','C2','C5','C6','C7'),]
C7 = genesets[genesets$gs_cat %in% c('C7'),]
C6 = genesets[genesets$gs_cat %in% c('C6'),]
C5 = genesets[genesets$gs_cat %in% c('C5'),]
C2 = genesets[genesets$gs_cat %in% c('C2'),]
H = genesets[genesets$gs_cat %in% c('H'),]


####直接循环每个
library(ggplot2)
library(ggpubr)
library(Seurat)
#BiocManager::install('limma')
library(limma)
a = list(H=H,C2=C2,C5=C5,C6=C6,C7=C7)
geneset_name=a[[1]]
Idents(sce) = sce$type
###需要根据目的来修改用于GSVA的对象
run_gsva_analysis <- function(geneset_name) {
  genesets = geneset_name
  genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
  genesets <- split(genesets$gene_symbol, genesets$gs_name)
  
  gsva_data <- as.data.frame(sce@assays$RNA$data)
  result <- gsvaParam(as.matrix(gsva_data), 
                 genesets,
                 kcdf="Gaussian")
  colnames(sce@meta.data)
  table(sce@meta.data$type)
  ###这里需要修改
  cluster <- sce@meta.data %>% dplyr::filter(type %in% c('MOS','Tissue')) %>% dplyr::arrange(type)
  
  ## limma 差异分析
  #需要运行一段时间
  use_gsva <- gsva(result)
  dim(use_gsva)
  head(use_gsva)[1:2,1:4]
  phe = sce@meta.data
  sce$type
  table(phe$type)
  p = identical(colnames(use_gsva),rownames(cluster));p
  if(!p) use_gsva = use_gsva[,match(rownames(cluster),colnames(use_gsva))]
  group <- c(rep("MOS", 1141), rep("Tissue", 1015)) %>% as.factor()
  desigN <- model.matrix(~ 0 + group) 
  colnames(desigN) <- levels(group)
  fit = lmFit(use_gsva, desigN)
  ##注意不要比反了，这里代表MOS（比较组）对比Tissue（对照）
  cont.matrix <- makeContrasts(contrasts = c('MOS-Tissue'), levels = desigN)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  diff <- topTable(fit2,adjust='fdr', coef=1, number=Inf)
  cluster2_diff <- na.omit(diff)
  cluster2_diff$pathway <- rownames(cluster2_diff)
  
  
  ### 筛选显著通路
  sig_cluster2_diff <- cluster2_diff %>% dplyr::filter(abs(logFC) > 0.2 & adj.P.Val < 0.05)
  
  draw_result <- sig_cluster2_diff
  library(dplyr)
  draw_result <- draw_result %>% dplyr::mutate(label = if_else(logFC > 0,"up","down")) %>% dplyr::arrange(logFC)
  
  draw_result$label <- factor(draw_result$label)
  p <- ggplot(draw_result,
              aes(x =logFC, y = pathway, fill = label)) + 
    geom_col() + 
    theme_bw() +
    theme(
      legend.position = 'none',
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = 'grey60',size = 1.1),
      axis.text = element_text(size = 12)
    )
  p
  
  #指定因子；调整顺序：
  draw_result$pathway <- factor(draw_result$pathway,levels = rev(draw_result$pathway))
  
  #先根据上下调标签拆分数据框：
  up <- draw_result[draw_result$label == 'up',]
  down <- draw_result[draw_result$label == 'down',]
  #添加上调pathway标签：
  p2 <- p +
    geom_text(data = up,
              aes(x = -0.05, y = pathway, label = pathway),
              size = 3.5,
              hjust = 1)+ #标签右对齐
    geom_text(data = down,
              aes(x = 0.05, y = pathway, label = pathway),
              size = 3.5,
              hjust = 0) #标签左对齐
  p2
  
  
  #install.packages('cols4all')
  library(cols4all) 
  c4a_gui()
  mycol <- c4a('bright',2)
  mycol
  p3 <- p2 +
    scale_x_continuous(breaks=seq(-4, 6, 2)) + #x轴刻度修改
    labs(x = 'Normalized Enrichment Score', y = ' ', title = ' ') + #修改x/y轴标签、标题添加
    theme(plot.title = element_text(hjust = 0.5, size = 14))+ #主标题居中、字号调整
    scale_fill_manual(values = mycol)
  p3
  
  ###个性化调节
  p4 <- p3 +
    geom_segment(aes(x = 0.5, y = 4, xend = 1 , yend = 4),  # 绘制向右的箭头
                 arrow = arrow(length = unit(0.03, "npc")),      # 箭头长度和类型
                 color = '#EE6677') + # 箭头颜色
    geom_text(aes(x = 0.9, y = 4 - 0.3, label = "Up in MOS", hjust = 1), # 添加注释
              color = '#EE6677', size = 5) + # 注释颜色和大小
    geom_segment(aes(x = -0.5, y = 4, xend = -1, yend = 4 ), # 绘制向下的箭头
                 arrow = arrow(length = unit(0.03, "npc")), # 箭头长度和类型
                 color = '#4477AA') + # 箭头颜色
    geom_text(aes(x = -0.5, y = 4 - 0.3, label = "Up in Tissue", hjust = 1), # 添加注释
              color = '#4477AA', size = 5) # 注释颜色和大小
  
  p4
  
  return(list(plot=p4, table=sig_cluster2_diff))
}


# results_list <- lapply(a, run_gsva_analysis)
# H_table = results_list$H$table
# H_plot = results_list$H$plot
# H_plot

results_list = run_gsva_analysis(a[[1]])
H_table = results_list$table
H_plot = results_list$plot
H_plot


rownames(H_table)
Hallmark = genesets[genesets$gs_name %in% rownames(H_table),]
table(Hallmark$gs_name)


###热图
library(pheatmap)
library(stringr)
data <- use_gsva[match(rownames(sig_cluster2_diff),rownames(use_gsva)),]

#对数据进行排序
identical(rownames(cluster),colnames(data))
Group <- factor(cluster$type,levels = c("MOS","Tissue"))
sort_data <- order(Group)
n <- data[,sort_data]
cluster <- cluster[sort_data,]

# 调整颜色梯度
range(data)
breaksList = seq(-0.5, 0.5, by = 0.1)
colors <- colorRampPalette(c("#336699", "white", "tomato"))(length(breaksList))

#创建列和行注释
annCol <- data.frame(group = cluster$type,
                     row.names = colnames(data),
                     stringsAsFactors = FALSE)
#annRow <- data.frame(row.names = rownames(data))
pheatmap(data,
         annotation_col = annCol,
         #annotation_row = annRow,
         color = colors,
         breaks = breaksList,
         cluster_rows = T,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         #gaps_col = cumsum(table(annCol$Type)),  # 使用排序后的列分割点
         #gaps_row = cumsum(table(annRow$Methods)), # 行分割
         fontsize_row = 6,
         fontsize_col = 6,
         annotation_names_row = FALSE
)





###GSVA分析
#看看选择的是：①基于单细胞表达量矩阵的gsea和gsva（需要生物学功能数据库）
#②基于细胞亚型表达量矩阵的gsea和gsva（需要生物学功能数据库）

####这里是第②种####
#多个分组后的gsva值，很难limma进一步得到logFC和p值，并不是两两分组；如果是两两分组，完全可以选择第一种方式
#按照细胞亚群平均基因
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(limma)
#BiocManager::install('clusterProfiler')
#genesets <- msigdbr(species = "Homo sapiens") 
#genesets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
genesets <- msigdbr(species = "Homo sapiens",category = 'H') 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)
?AverageExpression
Idents(sce) = sce$celltype
expr <- AverageExpression(sce, assays = "RNA", layer = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #选取非零基因
expr <- as.matrix(expr)
head(expr)

# gsva默认开启全部线程计算
?gsvaParam
gsva.res <- gsvaParam(expr, genesets) 
saveRDS(gsva.res, "gsva.res.rds")
#gsva.res = readRDS( "gsva.res.rds")
expr_geneset <- gsva(gsva.res)
dim(expr_geneset)
head(expr_geneset)[1:2,1:4]
write.csv(expr_geneset, "gsva_res.csv", row.names = F)
gsva_d = expr_geneset[sample(nrow(expr_geneset),30),]
pheatmap::pheatmap(gsva_d, show_colnames = T, 
                   scale = "row",angle_col = "45",
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50))


# 气泡图
library(reshape2)
gsva_long <- melt(gsva_d, id.vars = "Genesets")

# 创建气泡图
ggplot(gsva_long, aes(x = Var2, y = Var1, size = value, color = value)) +
  geom_point(alpha = 0.7) +  # 使用散点图层绘制气泡，alpha设置点的透明度
  scale_size_continuous(range = c(1, 6)) +  # 设置气泡大小的范围
  theme_minimal() + 
  scale_color_gradient(low='#427183',high='#D2D470') +
  labs(x = "", y = "Geneset", size = "GSVA Score")+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))

