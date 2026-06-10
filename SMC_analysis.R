##### Convert snRNAseq data file from python "anndata" to R "seurat" #####
library(sceasy)
library(LoomExperiment)
library(anndata)
library(reticulate)
library(Seurat)
library(tidyverse)
library(cowplot)

py_require("anndata") # works as "library" for python so the necessary packages are installed
sceasy::convertFormat("data/Chow_HFD.h5ad", from="anndata", to="seurat",
                      outFile="data/Chow_HFD.rds")
snData <- readRDS(file = "data/Chow_HFD.rds")

# As counts are not normalised in the converted object --> normalize first. As there seems to be an issue with the Seurat package and the normalize function, a fixed version is installed
devtools::install_github("satijalab/seurat", ref = "fix/v.5.3.1")
snData <- NormalizeData(snData)

##### Re-map celltypes based on the previous analysis of cell types #####

# Select resolution of clusters
Idents(snData) <- snData@meta.data$leiden_res0_25
snData@meta.data$new_celltypes <- snData@meta.data$leiden_res0_25

# Correct cluster labelling of the Chow clusters
levels(snData@meta.data$new_celltypes[snData@meta.data$batch == "Chow"]) <- c("0", "1", "2", "4", "6", "3", "5", "7")
Idents(snData) <- snData@meta.data$new_celltypes
# In the Chow cluster the ACs are in one cluster with the ECs. Use another resolution to subcluster (map cluster 4 res0_5 as cluster 5 on new_celltypes)
snData@meta.data$new_celltypes[snData@meta.data$batch == "Chow" & snData@meta.data$leiden_res0_5 == "4"] <- "5"
# In the HFD_L1 cluster the FBMC cluster (7) is merged with cluster 4. Use another resolution to subcluster (map cluster 7 res0_5 as cluster 7 on new_celltypes)
snData@meta.data$new_celltypes[snData@meta.data$batch == "HFD_L1" & snData@meta.data$leiden_res0_5 == "7"] <- "7"

Idents(snData) <- snData@meta.data$new_celltypes

# Rename Clusters according to previous cell type identification
snData <- RenameIdents(object = snData,
                       "0"   = "SMC",
                       "1"   = "SMC",
                       "2"   = "EC",
                       "3"   = "Immune",
                       "4"   = "FB",
                       "5"   = "AC",
                       "6"   = "SMC",
                       "7"   = "FBMC")
snData@meta.data$new_celltypes <- Idents(snData)

# Reorder identities
levels(x = snData) <- c("SMC", "FBMC", "FB", "EC", "AC", "Immune")

snData@meta.data$new_celltypes <- Idents(snData)

##### Find differences in expression due to treatment (separately for SMC and FBMC clusters)
# Combine celltype and treatment label
snData$new_celltypes.batch <- paste(snData$new_celltypes, snData$batch, sep = "_")

# Find markers
Idents(snData) <- "new_celltypes.batch"
SMC.diff <- FindMarkers(snData, ident.1 = "SMC_HFD", ident.2 = "SMC_HFD_L1")
FBMC.diff <- FindMarkers(snData, ident.1 = "FBMC_HFD", ident.2 = "FBMC_HFD_L1")

# Add gene symbols to the DE table
annotations <- read.csv("data/annotation.csv")
SMC.diff <- SMC.diff %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

FBMC.diff <- FBMC.diff %>%
  rownames_to_column(var = "gene") %>%
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

# Reorder columns and sort by padj
SMC.diff <- SMC.diff[, c(1, 3:5, 2, 6:7)]
SMC.diff <- SMC.diff %>%
  dplyr::arrange(p_val_adj)

FBMC.diff <- FBMC.diff[, c(1, 3:5, 2, 6:7)]
FBMC.diff <- FBMC.diff %>%
  dplyr::arrange(p_val_adj)

write.csv2(x = SMC.diff, file = "./results/SMC-HFDvsL1.csv")
write.csv2(x = FBMC.diff, file = "./results/FBMC-HFDvsL1.csv")

# >>> Filter list of DEGs in MS Excel and plot selected genes here afterwards <<< #

#Plot selected significantly DEGs by treatment for SMC and FBMC as ridgeplots
subNoChow <- subset(snData, subset = batch == "Chow", invert = T)
# SMC contractile
RidgePlot(object = subNoChow, features = c("Myh11", "Acta2", "Tagln", "Myl9", "Cnn1"), idents = "SMC", group.by = "batch")
# SMC ECM
RidgePlot(object = subNoChow, features = c("Col4a3", "Col15a1", "Col4a6", "Col1a2", "Fn1", "Col3a1", "Col1a1", "Bgn"), idents = "SMC", group.by = "batch")
# FBMC contractile
RidgePlot(object = subNoChow, features = c("Myh11", "Acta2", "Tagln"), idents = "FBMC", group.by = "batch")
# FBMC ECM
RidgePlot(object = subNoChow, features = c("Fn1","Bgn", "Col3a1"), idents = "FBMC", group.by = "batch")