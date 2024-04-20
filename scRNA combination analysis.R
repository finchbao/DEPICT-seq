
# -------------------------------------------------------------------------
# important package
# -------------------------------------------------------------------------

library(Seurat)
library(ggplot2)
library(ggsci)
library(dplyr)
library(reshape2)
library(stringi)
library(ggsignif)

# -------------------------------------------------------------------------
# Read PBMC 1K
# -------------------------------------------------------------------------

pure_pbmc = Read10X('~/pbmc/')
pbmc <- CreateSeuratObject(counts = pure_pbmc,min.cells = 5, min.features = 500)
pbmc$group = 'pbmc'

# -------------------------------------------------------------------------
# Read 10X Jurkat 
# -------------------------------------------------------------------------

Jurkat_data = Read10X(data.dir = "~/Jurkat/")
Jurkat_data_matrix = as.matrix(Jurkat_data)
set.seed(43210) # set random seed
Jurkat_data_matrix_sample = Jurkat_data_matrix[,sample(1:3258,size = 50)] # sample 50 Jurkat cells
g2s = fread('/data/kxbao/20240316Library_analysis/g2s_jurkat_gencode.txt')
colnames(g2s) <- c("geneid","symbol")

# -------------------------------------------------------------------------
# Read Depict-seq mtx
# -------------------------------------------------------------------------

raw_mtx = read.table('~/depict_mtx.txt',comment.char = '#',header=T,row.names = 1) %>% .[,-c(1:4)]
g2s_my = fread('~/g2s_gencode.txt')
colnames(g2s_my) <- c("geneid","symbol")
raw_mtx$geneid = row.names(raw_mtx)

depict_mtx = merge(g2s_my,raw_mtx,by.x = 'geneid',by.y = 'geneid',all = T)
depict_mtx =depict_mtx[,-1]
depict_mtx_sum <- aggregate(.~symbol,FUN=sum,depict_mtx) # combine the reads from same gene
row.names(my_mtx_sum) = my_mtx_sum[,1]
my_mtx_sum =my_mtx_sum[,-1]
saveRDS(my_mtx_sum,'~/depict_dupensefliter_mtx.rds')


# -------------------------------------------------------------------------
# Combine Analysis
# -------------------------------------------------------------------------

# Combine Jurkat from DEPICT-seq and 10X

Jurkat_data_matrix_sample = as.data.frame(Jurkat_data_matrix_sample)
Jurkat_matrix = merge(Jurkat_data_matrix_sample,depict_mtx_sum,by = 0,all=T)
row.names(Total_matrix) = Total_matrix[,1]
Total_matrix = Total_matrix[,-1]
dataan <-  as(as.matrix(Total_matrix), "dgCMatrix")
Jurkat_10X <- CreateSeuratObject(counts = dataan)
Jurkat_10X$group =c(rep('Jurkat',50),rep('DEPICT',9),rep('Live',9))

# Combine Jurkat to PBMC 1K

combined = merge(pbmc,Jurkat_10X)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined,npcs = 5)
combined <- FindNeighbors(combined, dims = 1:5, reduction = "pca")
combined <- FindClusters(combined, resolution = 0.5, cluster.name = "unintegrated_clusters")
combined <- RunTSNE(combined, reduction = "pca", dims = 1:5, reduction.name = "tsne.unintegrated")
DimPlot(combined,reduction = "tsne.unintegrated",group.by = c("group", "seurat_clusters"),label = T)

# -------------------------------------------------------------------------
# Marker
# -------------------------------------------------------------------------
marker = c("LYZ", # Mono
           "CD14", # Mono
           "FCGR3A", # Mono
           "CD3E", # T
           "CD8A", # CD8 T
           'CCR7'# Naive T
           'MZB1', #Jurkat
           "GNLY", # NK
           'KLRF1', # NK
           "MS4A1", # B
)

FeaturePlot(combined, reduction = "tsne.unintegrated",features = marker)

FeaturePlot_scCustom(seurat_object = combined, features = marker, num_columns = 5)

combined <- RenameIdents(combined, 
                         `0` = "CD14+ Mono",
                         `1` = "Memory CD4 T", 
                         `2` = "B Cell",
                         `3` = "Naive CD4 T", 
                         `4` = "CD8 T", 
                         `5` = "Jurkat", 
                         `6` = "FCGR3A+ Mono",
                         `7` = "NK")

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)


VlnPlot(combined, features = marker,
        col = my36colors,
        stack = TRUE, 
        sort = TRUE, 
        flip = TRUE) +
  theme(legend.position = "none")

DimPlot(combined, reduction = "tsne.unintegrated", label = F)+  scale_color_frontiers()+theme_classic()
