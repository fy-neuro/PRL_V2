library(Seurat)
library(CellChat)
library(patchwork)
library(gridExtra)
library(writexl)

setwd("/storage/liuxiaodongLab/liaozizhuo/Projects/fanyi/code/")

job_id <- "2"

sobj <- readRDS("/storage/liuxiaodongLab/fanyi/FY/PRL_V2/20250123-24-Annotation-third-MNN/Documents/sobj_3rd_cleaned.rds")
sobj <- subset(sobj,subset=subclass != "T cells" & subclass != "BAMs" & subclass != "Vascular")
unique(sobj$Stage)
sobj <- subset(sobj, Stage %in% "Control")
unique(sobj$Stage)
Idents(sobj) <- "celltype"

data.input <- GetAssayData(sobj, slot = "data", assay = "RNA")
md <- sobj@meta.data
meta <- md[ , c("Sample", "subclass")]
colnames(meta) <- c("samples", "labels")
meta$samples <- factor(meta$samples)
# meta$labels <- paste0("c", meta$labels)
meta$labels <- factor(meta$labels)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "RNA")
cellchat

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 1) 
# Increase the maximum allowed size of global variables to 1.5 GB
options(future.globals.maxSize = 3 * 1024^3)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = 100)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

# Define cluster names and colors
cluster_names <- c("ID2", "PV", "SST", "VIP", 
                   "L2/3-CUX2", "L3-5-RORB", "L5/6-THEMIS", "L5/6-TLE4", 
                   "Astro-PLCG1", "Astro-SERPINI2", "Microglia", 
                   "Microglia-IL1B", "T cells", "BAMs", "Oligos", "OPCs", "OPCs-MBP", "Vascular")
celltype_colors <- c("#65b38a", "#F3746C", "#ec3f4c", "#a6cee3", "#588198", 
                     "#FA8737", "#E9AA44", "#C45338", "#17DAC6", "#BB4A94", 
                     "#4993fa", "#e5f5e0", "#a1d99b", "#ffcc00", "#8c564b", 
                     "#c49c94", "#9467bd", "#d62728")

# Set cluster names in the CellChat object
levels(cellchat@idents) <- cluster_names
names(celltype_colors) <- cluster_names
vertex_colors <- celltype_colors[levels(cellchat@idents)]

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, color.use = vertex_colors, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, color.use = vertex_colors, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Open a PDF device
pdf(paste0("../plots/netVisual_circle_plot_Control_", job_id, ".pdf"))

# Create the plot
netVisual_circle(cellchat@net$count, 
                 vertex.weight = rowSums(cellchat@net$count), 
                 weight.scale = TRUE, 
                 color.use = vertex_colors,
                 label.edge = FALSE, 
                 title.name = "Number of interactions")

# Close the PDF device
dev.off()

cellchat@netP$pathways

# List of pathways
pathways <- cellchat@netP$pathways

# Create a PDF to store all plots
pdf(paste0("../plots/all_pathways_plots_Control_", job_id, ".pdf"))

# Loop through each pathway and create a plot
for (pathway in pathways) {
  pathways.show <- c(pathway)
  
  # Generate the plot with the title corresponding to the pathway
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", title = pathway, color.use = vertex_colors)
  
  # Add the title manually
  title(main = pathway, cex.main = 1.5)  # You can adjust cex.main for title size
}

# Close the PDF
dev.off()

df.net <- subsetCommunication(cellchat)
write_xlsx(df.net, paste0("../output/Ligand-Receptor_interaction_Control_", job_id, ".xlsx"))

# Compute network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat)

# Generate signaling role heatmap
pdf(paste0("../plots/netAnalysis_signalingRole_heatmap_Control_", job_id, ".pdf"))
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", color.use = vertex_colors, font.size = 5)
dev.off()
