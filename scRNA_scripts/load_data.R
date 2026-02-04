####单个10X样本读取####
#获取目录中的所有文件名
file.name <- list.files("GSE239676_RAW/")
#查找特定类型的文件名
barcode <- file.name[grep("barcodes.tsv.gz$", file.name)]
feature <- file.name[grep("features.tsv.gz$", file.name)]
matrix <- file.name[grep("matrix.mtx.gz$", file.name)]
#读取并加载文件内容
data1 <- as.data.frame(readMM(paste0("GSE239676_RAW/", matrix[1])))
data1[1:6,1:6]
data2 <- fread(paste0("GSE239676_RAW/", barcode[1]), header = FALSE)
data2[1:6,]
data3 <- fread(paste0("GSE239676_RAW/", feature[1]), header = FALSE)
data3[1:6,]


colnames(data1) <- data2$V1 #细胞barcode
data1$gene_id <- data3$V2 #基因名
data1 <- data1[!duplicated(data1$gene_id),] # 去除重复基因
rownames(data1) <- data1$gene_id
data1 <- data1 %>% dplyr::select(-gene_id) #把基因ID放到第一列
#使用原始（非标准化数据）初始化Seurat对象
sce <- CreateSeuratObject(counts = data1,   min.cells = 5,
                          min.features = 300)


####多个10X样本读取####
# 解压缩获取数据
# 一般下载下来的都是tar结尾的压缩文件
untar("GSE205013_RAW.tar",exdir = "GSE205013_RAW")
list.files("GSE205013_RAW/")
dir='GSE205013_RAW/' 
fs=list.files('GSE205013_RAW/','^GSM')
fs
library(tidyverse)
samples=str_split(fs,'_',simplify = T)[,1]

##处理数据，将原始文件分别整理为barcodes.tsv.gz，features.tsv.gz和matrix.mtx.gz到各自的文件夹
#批量将文件名改为 Read10X()函数能够识别的名字
if(F){
  lapply(unique(samples),function(x){
    # x = unique(samples)[1]
    y=fs[grepl(x,fs)]
    folder=paste0("GSE205013_RAW/", paste(str_split(y[1],'_',simplify = T)[,1:2], collapse = "_"))
    dir.create(folder,recursive = T)
    #为每个样本创建子文件夹
    file.rename(paste0("GSE205013_RAW/",y[1]),file.path(folder,"barcodes.tsv.gz"))
    #重命名文件，并移动到相应的子文件夹里
    file.rename(paste0("GSE205013_RAW/",y[2]),file.path(folder,"features.tsv.gz"))
    file.rename(paste0("GSE205013_RAW/",y[3]),file.path(folder,"matrix.mtx.gz"))
  })
}

samples=list.files( dir )
samples 
sceList = lapply(samples,function(pro){ 
  # pro=samples[1] 
  print(pro)  
  tmp = Read10X(file.path(dir,pro )) 
  if(length(tmp)==2){
    ct = tmp[[1]] 
  }else{ct = tmp}
  sce =CreateSeuratObject(counts =  ct ,
                          project =  pro  ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)
}) 

do.call(rbind,lapply(sceList, dim)) 

sce.all=merge(x=sceList[[1]],
              y=sceList[ -1 ],   add.cell.ids = samples ) 
names(sce.all@assays$RNA@layers)
sce.all[["RNA"]]$counts 
# Alternate accessor function with the same result
LayerData(sce.all,assay = "RNA", layer = "counts")
sce.all <- JoinLayers(sce.all)
dim(sce.all[["RNA"]]$counts )

as.data.frame(sce.all@assays$RNA$counts[1:10, 1:2])
head(sce.all@meta.data, 10)
table(sce.all$orig.ident) 


####单个h5样本读取####
pro = "train"
list.files("input/")
# [1] "GSM8128607_P1_B_filtered_feature_bc_matrix.h5"
sce=CreateSeuratObject(counts =  Read10X_h5("input/GSM8128607_P1_B_filtered_feature_bc_matrix.h5") ,
                       project =  pro ,
                       min.cells = 5,
                       min.features = 300,)
dim(sce)
#[1] 23062 11442


####多个h5数据样本读取####
#install.packages('hdf5r')
library(hdf5r)
library(tidyverse)
dir='GSE184242_RAW' 
samples=list.files('GSE184242_RAW/','^GSM')
samples

sceList = lapply(samples,function(pro){ 
  print(pro) 
  sce =CreateSeuratObject(counts = Read10X_h5( file.path(dir,pro)) ,
                          project =  gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',pro) ) ,
                          min.cells = 5,
                          min.features = 300 )
  return(sce)
})
sceList
samples

#整合数据
sce.all=merge(x=sceList[[1]],y=sceList[ -1 ],
              add.cell.ids =  gsub('_filtered_feature_bc_matrix.h5','',gsub('^GSM[0-9]*_','',samples)))
dim(sce.all) 
# [1] 27747 78306



####单个txt/csv/tsv数据读取###
pro = "train"
list.files("input/")
#[1] "GSE181919_UMI_counts.txt.gz"
ct=fread("input/GSE181919_UMI_counts.txt.gz",data.table = F)
rownames(ct) <- ct$V1
ct <- ct[,-1]
head(ct)[1:5,1:2]
#               AAACGGGCATGACGGA.1 AAAGATGAGCAGACTG.1
# RP11-34P13.7                   0                  0
# FO538757.2                     1                  3
# AP006222.2                     0                  0
# RP4-669L17.10                  0                  0
# RP11-206L10.9                  0                  0
sce=CreateSeuratObject(counts =  ct ,
                       project =  pro ,
                       min.cells = 5,
                       min.features = 300,)
dim(sce)


####多个txt/csv/tsv文件读取####
untar("GSE167297_RAW.tar",exdir = "GSE167297_RAW")
dir='GSE167297_RAW/' 
samples=list.files( dir )
#samples=list.files('GSE167297_RAW/','^GSM')
samples
# [1] "GSM5101019_Pt3_Superficial_CountMatrix.txt.gz" "GSM5101020_Pt3_Deep_CountMatrix.txt.gz"       
# [3] "GSM5101021_Pt4_Normal_CountMatrix.txt.gz" 

sceList = lapply(samples,function(pro){ 
  print(pro)
  ct=fread(file.path( dir ,pro),data.table = F)
  ct[1:4,1:4]
  rownames(ct)=ct[,1]
  ct=ct[,-1]
  sce=CreateSeuratObject(counts =  ct ,
                         project = gsub('_CountMatrix.txt.gz','',gsub('^GSM[0-9]*_','',pro) ) ,
                         min.cells = 5,
                         min.features = 300)
  
  return(sce)
})
sceList
samples

#整合数据
sce.all=merge(x=sceList[[1]],y=sceList[ -1 ],
              add.cell.ids =  gsub('_CountMatrix.txt.gz','',gsub('^GSM[0-9]*_','',samples)))
dim(sce.all)
# [1] 16806  7827



