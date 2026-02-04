
#BiocManager::install ('ReactomePA')
com_go_kegg_ReactomePA_human <- function(symbols_list ,pro){ 
  library(clusterProfiler)
  library(org.Hs.eg.db)
  #library(ReactomePA)
  library(ggplot2)
  library(stringr)
  
  # 首先全部的symbol 需要转为 entrezID
  gcSample = lapply(symbols_list, function(y){ 
    y=as.character(na.omit(AnnotationDbi::select(org.Hs.eg.db,
                                  keys = y,
                                  columns = 'ENTREZID',
                                  keytype = 'SYMBOL')[,2])
    )
    y
  })
  gcSample
  
  # 第1个注释是 KEGG 
  xx <- compareCluster(gcSample, fun="enrichKEGG",
                       organism="hsa", pvalueCutoff=0.05)
  dotplot(xx)  + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro,'_kegg.pdf'),width = 10,height = 8)
  
  # 第2个注释是 ReactomePA 
  xx <- compareCluster(gcSample, fun="enrichPathway",
                       organism = "human",
                  pvalueCutoff=0.05)
  dotplot(xx)  + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=50)) 
  ggsave(paste0(pro,'_ReactomePA.pdf'),width = 10,height = 8)
  
  # 然后是GO数据库的BP,CC,MF的独立注释
  # Run full GO enrichment test for BP 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont		   = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_BP_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_BP_cluster_simplified.csv'))
  # Run full GO enrichment test for CC 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont		   = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_CC_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_CC_cluster_simplified.csv'))
  
  # Run full GO enrichment test for MF 
  formula_res <- compareCluster(
    gcSample, 
    fun="enrichGO", 
    OrgDb="org.Hs.eg.db",
    ont		   = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # Run GO enrichment test and merge terms 
  # that are close to each other to remove result redundancy
  lineage1_ego <- simplify(
    formula_res, 
    cutoff=0.5, 
    by="p.adjust", 
    select_fun=min
  ) 
  pdf(paste0(pro,'_GO_MF_cluster_simplified.pdf') ,width = 15,height = 8)
  print(dotplot(lineage1_ego, showCategory=5) + theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
          scale_y_discrete(labels=function(x) str_wrap(x, width=50)) )
  dev.off() 
  write.csv(lineage1_ego@compareClusterResult, 
            file=paste0(pro,'_GO_MF_cluster_simplified.csv'))
  
}
