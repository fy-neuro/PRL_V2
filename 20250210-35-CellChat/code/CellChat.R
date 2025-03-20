library(dplyr)
library(Seurat)
library(patchwork)
library(readr)
library(cols4all)
library(ggplot2)
library(openxlsx)
library(tibble)
library(loupeR)
library(SeuratWrappers)
library(CellChat)


##Set path
setwd("/storage/liuxiaodongLab/fanyi/FY/PRL_V2/20250210-35-CellChat")
plot_path <- "./plot"
Document_path <- "./Documents"


###load data
sobj<-readRDS("/storage/liuxiaodongLab/fanyi/FY/PRL_V2/20250123-24-Annotation-third-MNN/Documents/sobj_3rd_cleaned.rds")

sobj$subcluster<- sobj$subcluster %>% as.character %>% as.factor
sobj$samples <- sobj$orig.ident %>% as.factor


# 创建函数简化重复操作
create_cellchat <- function(seurat_obj, celltype_column) {
  # 提取表达矩阵和metadata
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  meta <- seurat_obj@meta.data[, c("orig.ident", "Group", celltype_column,"samples")]
  colnames(meta)[3] <- "celltype"
  
  # 创建CellChat对象
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")
  
  # 添加细胞注释
  cellchat <- setIdent(cellchat, ident.use = "celltype")
  
  # 预处理
  cellchat@DB <- CellChatDB.human  # 根据实际物种选择
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  return(cellchat)
}


# 对于大数据集可启用并行
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 16 * 1024^3)


# 对每个注释层级进行分析
for (ct_level in c("class", "subclass", "subcluster")) {
  cat(ct_level,"\n")
  # 创建CellChat对象
  cellchat <- create_cellchat(sobj, ct_level)
  
  # 计算通讯概率
  cellchat <- computeCommunProb(cellchat)
  cellchat <- computeCommunProbPathway(cellchat)
  
  saveRDS(cellchat,paste0("./Documents/Cellchat_Object_",ct_level,".rds"))
  # 可视化
  pathways.show <- cellchat@netP$pathways
  p1 <- netVisual_heatmap(cellchat)
  p2 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize)
  
  # 保存结果
  pdf(paste0("../plot/CellChat_", ct_level, ".pdf"), width = 12, height = 8)
  print(p1)
  print(p2)
  dev.off()
}


# 拆分Control和PD组
control_cells <- Cells(sobj)[sobj$Group == "Control"]
pd_cells <- Cells(sobj)[sobj$Group == "PD"]

# 创建分组对象
cellchat_control <- create_cellchat(sobj[, control_cells], "subcluster")
cellchat_pd <- create_cellchat(sobj[, pd_cells], "subcluster")

# 计算通讯概率（需要分别计算）
cellchat_control <- computeCommunProb(cellchat_control)
cellchat_pd <- computeCommunProb(cellchat_pd)

# 合并对象进行比较
cellchat_list <- list(Control = cellchat_control, PD = cellchat_pd)
saveRDS(cellchat_list,paste0("./Documents/Cellchat_Object_",ct_level,".rds"))

cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

# 比较总通讯强度
pdf("../plot/Group_comparison.pdf", width = 15, height = 10)
print(compareInteractions(cellchat, measure = "count"))
print(compareInteractions(cellchat, measure = "weight"))
dev.off()


# 前20通路分析（以PD组为例）
cellchat_pd <- computeNetSimilarity(cellchat_pd)
cellchat_pd <- netAnalysis_computeCentrality(cellchat_pd)

top20_pathways <- rankNet(cellchat_pd, mode = "comparison")$pathways[1:20]
pdf("./plot/Top20_Pathways_PD.pdf", width = 12, height = 8)
for (pathway in top20_pathways) {
  print(netVisual_heatmap(cellchat_pd, signaling = pathway))
}
dev.off()

# 前50通路差异比较
pdf("./plot/Top50_Differential_Pathways.pdf", width = 18, height = 12)
print(netVisual_heatmap(cellchat, comparison = c(1,2), 
                        signaling = cellchat@netP$pathways[1:50],
                        color.heatmap = "RdBu"))
dev.off()

