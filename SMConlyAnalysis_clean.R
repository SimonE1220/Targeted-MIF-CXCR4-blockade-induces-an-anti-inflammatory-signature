# Subclustering of SMC/FBMC clusters
# snData refers to the fully analyzed dataset

library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(ggplot2)
library(viridis) # necessary to use the viridis colour scheme
library(mascarade) # Allows to highlight clusters by outlines
library(purrr)
library(furrr)
library(plotly)


# Create SMC/FBMC subset and prepare it for analysis ----------------------

# Create subset with only SMC and FBMC cluster
SMC.only <- subset(snData, idents = c("SMC", "FBMC"))

# Briefly check whether cell cycle or mitochondrial genes have a relevant influence on PCA analysis (check for full data set showed no relevant influence)

# Load cell cycle markers to assess effects of cell cycle
load("data/ccMarkers_MusMusculus.RData") # Creates vectors "g2m_genes" and "s_genes"
SMC.only <- CellCycleScoring(SMC.only,
                             g2m.features = g2m_genes,
                             s.features = s_genes)

#Calculate percentage of mitochondrial genes to check whether it should be regressed out
SMC.only <- PercentageFeatureSet(SMC.only, pattern = "^mt-", col.name = 'percent.mt')

# Assess the effect of mitochondrial genes on clustering (mitoRatio)
# Check quartile values
summary(SMC.only@meta.data$percent.mt)

# Turn mitoRatio into categorical factor vector based on quartile values
SMC.only@meta.data$mitoFr <- cut(SMC.only@meta.data$percent.mt,
                                 breaks = c(-Inf, 0.3276, 0.8148, 1.6582, Inf),
                                 labels = c("Low", "Medium", "Medium high", "High"))

# Perform a quick and dirty PCA for testing of mitochondrial RNA and cell cycle effects
SMC.only <- NormalizeData(SMC.only)
SMC.only <- FindVariableFeatures(SMC.only,
                                 selection.method = "vst",
                                 nfeatures = 2000,
                                 verbose = TRUE)
SMC.only <- ScaleData(SMC.only)
SMC.only <- RunPCA(SMC.only)

# Check potential influence of confounders on clustering
DimPlot(SMC.only,
        reduction = "pca",
        group.by = "Phase",
        split.by = "batch")

DimPlot(SMC.only,
        reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")

DimPlot(SMC.only,
        reduction = "pca",
        group.by = "mitoFr",
        split.by = "mitoFr")

DimPlot(SMC.only,
        reduction = "pca",
        group.by = "mitoFr",
        split.by = "batch")

# There is no notalbe effect of the cell cycle score on PCA but a clear effect of mitochondrial genes on PCA (no batch effect but the same effect within each batch)
# REGRESS OUT MITOCHONDRIAL GENES

# SCTransformation and Integration ----------------------------------------
# As the Seurat object is in the v3 format, it will be converted to the new format (v5 at the moment of writing)
SMC.only <- UpdateSeuratObject(object = SMC.only)

#Split object, SCTransform, and IntegrateLayers
SMC.only[["RNA"]] <- split(SMC.only[["RNA"]], f = SMC.only$batch)
SMC.only <- SCTransform(object = SMC.only, vars.to.regress = "percent.mt")
SMC.only <- RunPCA(object = SMC.only, npcs = 50, assay = "SCT")

SMC.only <- IntegrateLayers(object = SMC.only, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", normalization.method = "SCT")

# Cluster data
SMC.only <- FindNeighbors(object = SMC.only, reduction = "integrated.rpca", dims = 1:10)
SMC.only <- FindClusters(object = SMC.only, resolution = c(0.1, 0.2, 0.4, 0.8, 1.5), cluster.name = c("rpca_clusters_0.1", "rpca_clusters_0.2", "rpca_clusters_0.4", "rpca_clusters_0.8", "rpca_clusters_1.5"))

# Generate UMAPs to inspect overall appearance
SMC.only <- RunUMAP(object = SMC.only, reduction = "integrated.rpca",    dims = 1:10, reduction.name = "umap.rpca")

DefaultAssay(SMC.only) <- "RNA"
SMC.only <- JoinLayers(SMC.only)
Idents(SMC.only) <- SMC.only$newCellred


# Check which resolution to take --------------------------

# Use the "mascarade" package to show the "old" cluster outlines to find a resolution that covers the old FBMC cluster best.
# As default, idents are used for the outline
maskTable <- generateMaskSeurat(SMC.only, reduction = "umap.rpca")

DimPlot(SMC.only,
        reduction = "umap.rpca",
        label = T,
        label.size = 10,
        pt.size = 0.5,
        repel = TRUE,
        group.by = "rpca_clusters_0.8",
 ) +
   fancyMask(maskTable, ratio = 1)

# CONTINUE WITH RPCA 0.8
Idents(SMC.only) <- "rpca_clusters_0.8"

# Create module scores ----------------------------------------------------
# The module score needs normalized and scaled data.
DefaultAssay(SMC.only) <- "SCT"

# Use a reduced marker list for the module score calculation
# Combined scRNAseq markers
sc.Marker.short <- list(contrVSMC     = c("Myh11", "Des", "Acta2", "Tagln", "Cnn1", "Myocd"),
                        modVSMC       = c("Col1a1", "Fn1", "Itga5", "Rbp1", "Spp1", "S100a4"),
                        myofibVSMC    = c("Pdgfrb", "Fn1", "Col1a1", "Lum", "Dcn", "Bgn"),
                        mesenchVSMC   = c("Cd34", "Cd44", "Eng", "Ly6a", "Klf4"),
                        macroVSMC     = c("Cd68", "Itgam", "Adgre1", "Lamp2", "Lgals3"),
                        foamVSMC      = c("Serpine1", "Cfh", "Cd36", "Thy1"),
                        semVSMC       = c("Vcam1", "Ly6c1"),
                        osteoVSMC     = c("Runx2", "Sox9", "Alpl", "Bmp2", "Msx2", "Spp1"),
                        adipoVSMC     = c("Ucp1", "Cfd", "Adipoq", "Lep", "Pparg"),
                        Pericyte      = c("Pdgfrb", "Ndufa4l2", "Notch3", "Cspg4", "Rgs5", "Gucy1b1"),
                        Prol          = c("Mki67", "Pcna", "Cenpa", "Aurkb", "Mcm2", "Top2a"))

# Calculate module scores. This will create metadata columns for each module (cell type)
SMC.only <- AddModuleScore(
  object = SMC.only,
  features = sc.Marker.short,
  name = "cellSignature"
)

# Signature renaming ------------------------------------------------------
# Give signatures speaking names
SMC.only$sig.contrVSMC <- SMC.only$cellSignature1
SMC.only$cellSignature1 <- NULL
SMC.only$sig.modVSMC <- SMC.only$cellSignature2
SMC.only$cellSignature2 <- NULL
SMC.only$sig.myofibVSMC <- SMC.only$cellSignature3
SMC.only$cellSignature3 <- NULL
SMC.only$sig.mesenchVSMC <- SMC.only$cellSignature4
SMC.only$cellSignature4 <- NULL
SMC.only$sig.macroVSMC <- SMC.only$cellSignature5
SMC.only$cellSignature5 <- NULL
SMC.only$sig.foamVSMC <- SMC.only$cellSignature6
SMC.only$cellSignature6 <- NULL
SMC.only$sig.semVSMC <- SMC.only$cellSignature7
SMC.only$cellSignature7 <- NULL
SMC.only$sig.osteoVSMC <- SMC.only$cellSignature8
SMC.only$cellSignature8 <- NULL
SMC.only$sig.adipoVSMC <- SMC.only$cellSignature9
SMC.only$cellSignature9 <- NULL
SMC.only$sig.Pericyte <- SMC.only$cellSignature10
SMC.only$cellSignature10 <- NULL
SMC.only$sig.Prol <- SMC.only$cellSignature11
SMC.only$cellSignature11 <- NULL
SMC.only$sig.senescentVSMC <- NULL

# Module score and marker gene plotting ---------------------------------------------------
DefaultAssay(SMC.only) <- "RNA"

# Create list of cell type marker genes
sc.Marker <- list(contrVSMC     = c("Myh11", "Des", "Acta2", "Tagln", "Cnn1", "Myocd", "Smtn", "Cald1", "Itga8", "Lmod1", "Kcnmb1", "Kcna5", "Synm", "Mylk", "Myl9"),
                  modVSMC       = c("Col1a1", "Fn1", "Tnc", "Itga5", "Itga9", "Rbp1", "Spp1", "Lum", "S100a4", "Mmp9", "Lcn2", "Hmga1"),
                  myofibVSMC    = c("Pdgfrb", "S100a4", "Fn1", "Col1a1", "Lum", "Dcn", "Bgn", "Tcf21", "Tnfrsf11b", "Vim", "Pxn"),
                  mesenchVSMC   = c("Cd34", "Cd44", "Eng", "Ly6a", "Klf4"),
                  macroVSMC     = c("Cd68", "Itgam", "Adgre1", "Lamp2", "Ccl2", "Csf2ra", "Hmox1", "Thbs1", "Lgals3"),
                  foamVSMC      = c("Serpine1", "Cfh", "Cd36", "Thy1"),
                  semVSMC       = c("Vcam1", "Ly6c1"),
                  osteoVSMC     = c("Runx2", "Sox9", "Alpl", "Bmp2", "Msx2", "Sp7", "Bglap", "Col2a1", "Col10a1", "Spp1"),
                  adipoVSMC     = c("Ucp1", "Cfd", "Prdm16", "Ppargc1a", "Adipoq", "Tmem26", "Trpv1", "Cebpa", "Lep", "Pparg"),
                  senescentVSMC = c("Cdkn2a", "Cdkn1a"), # only one of the two markers is found so it makes no sense to create a module score for them --> remove in short list
                  Pericyte      = c("Pdgfrb", "Ndufa4l2", "Notch3", "Cspg4", "Des", "Rgs5", "Gucy1b1"),
                  Prol          = c("Mki67", "Pcna", "Cenpa", "Aurkb", "Mcm2", "Top2a"))

# Plot module scores for cell types
ct.signatures <- c("sig.contrVSMC",
                   "sig.modVSMC",
                   "sig.myofibVSMC",
                   "sig.mesenchVSMC",
                   "sig.macroVSMC",
                   "sig.foamVSMC",
                   "sig.semVSMC",
                   "sig.osteoVSMC",
                   "sig.adipoVSMC",
                   "sig.senescentVSMC",
                   "sig.Pericyte",
                   "sig.Prol")

FeaturePlot(object = SMC.only,
            reduction = "umap.rpca",
            features = ct.signatures,
            order = TRUE,
            max.cutoff = "q95",
            pt.size = 0.5,
            label = TRUE,
            repel = TRUE) &
  scale_colour_gradient2(      # Add this colour scale to also make negative score values visible
    low = "blue",
    mid = "lightgrey",
    high = "red",
    midpoint = 0
  )

DotPlot(SMC.only,
        features = ct.signatures
        ) +
  scale_color_viridis(direction = 1, option = "D") +
  theme(axis.text.x = element_text(angle = 90))

# Plot marker gene expression to check its correlation with the clusters (check cell type wise)
FeaturePlot(object = SMC.only,
            reduction = "umap.rpca",
            features = sc.Marker$contrVSMC,
            cols = c("lightgrey", "red"),
            order = TRUE,
            min.cutoff = "q10",
            max.cutoff = "q95",
            pt.size = 0.5,
            label = TRUE,
            repel = TRUE)

DotPlot(SMC.only,
        features = sc.Marker$contrVSMC
        ) +
  scale_color_viridis(direction = 1, option = "D") +
  theme(axis.text.x = element_text(angle = 90))

# Merge similar clusters and assign cell types --------------------------------------------------
SMC.only <- RenameIdents(SMC.only,
                         "4"  = "SMC1",
                         "10" = "SMC1",
                         "0"  = "SMC1",
                         "2"  = "SMC1",
                         "1"  = "SMC2",
                         "6"  = "SMC2",
                         "7"  = "SMC2",
                         "11" = "SMC2",
                         "12" = "SMC2",
                         "20" = "SMC2",
                         "13" = "SMC3",
                         "14" = "SMC3",
                         "15" = "SMC3",
                         "17" = "SMC3",
                         "18" = "SMC3",
                         "19" = "SMC3",
                         "22" = "SMC3",
                         "24" = "SMC3",
                         "5"  = "SMC4",
                         "23" = "SMC5",
                         "8"  = "SMC6",
                         "9"  = "SMC7",
                         "25" = "SMC8",
                         "21" = "SMC9",
                         "16" = "FBMC")

# DE analysis for individual clusters between conditions ------------------
# Add celltype-batch correlation to metadata
SMC.only$celltype <- Idents(SMC.only)
SMC.only$celltype_batch <- paste(SMC.only$celltype, SMC.only$batch, sep = "_")
Idents(SMC.only) <- "celltype_batch"
cluster_list <- levels(SMC.only$celltype)

# Generate a grid to allow to pass multiple arguments to the mapping function
parameter_grid <- expand.grid(
  cluster = cluster_list,
  cond_1  = "HFD_L1",
  cond_2  = "HFD",
  stringsAsFactors = FALSE
)

# 1. Set up parallel processing
plan(multisession, workers = 12)

# 2. Define the function: IMPORTANT the variable names have to match the columns of the grid function
get_DEGs <- function(cluster, cond_1, cond_2){
  message(paste("Starting parallel processing for Cluster", cluster, "between conditions", cond_1, "and", cond_2))
  
  group1 <- paste(cluster, cond_1, sep = "_")
  group2 <- paste(cluster, cond_2, sep = "_")
  
  markers <- FindMarkers(SMC.only,
                         ident.1 = group1,
                         ident.2 = group2,
                         test.use = "MAST"
                        )
  
  markers %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    mutate(cluster_id = cluster) %>%
    select(cluster_id, everything()) # Puts cluster_id to the front
}

HFD_L1.vs.HFD <- future_pmap_dfr(#pmap can handle multiple arguments
  .l = parameter_grid,
  .f = get_DEGs,
  .options = furrr_options(seed = TRUE)
)  # Runs the loop in parallel instead of sequential with "map_dfr"

# Return to sequential processing
plan(sequential)

# Plot frequencies to check how relative cluster size changes due to treatment --------------------------------------------------------
Idents(SMC.only) <- "celltype"
# Plot percentage of cells per cluster
# Extract required data
clusterInfo <- FetchData(SMC.only, vars = c("ident", "batch"))

# Check cluster sizes in a table
clusterInfo %>% dplyr::count(ident, batch) %>% tidyr::pivot_wider(names_from = ident, values_from = n)

# Figure exports ---------------------------------------------------
Idents(SMC.only) <- "celltype"

# Plot bars grouped by sample as percentage:
# 'group' defines how to group for determination of percentage
# 'geom_bar(aes(y = ..prop..), position = "dodge", stat = "count")': ..prop.. for proportion determination by stat = "count", "dodge" places the bars next to each other. Otherwise they are stacked on top
ggplot(clusterInfo,
       aes(x = ident,
           fill = batch,
           group = batch)) +
  geom_bar(aes(y = after_stat(prop)),
           position = "dodge",
           stat = "count") +
  scale_y_continuous(labels = scales::percent) +
  ylab("Percentage") +
  xlab("Cluster") +
  ggtitle("Distribution of SMCs/FBMCs between clusters") +
  theme(plot.title = element_text(hjust = 0.5))

# UMAP
DimPlot(SMC.only,
        reduction = "umap.rpca",
        label = F,
        label.size = 10,
        pt.size = 0.5,
        repel = TRUE,
        split.by = "batch"
) +
  labs(title = "UMAP of subclustered SMCs/FBMCs - split by condition", x = "UMAP 1", y = "UMAP 2") +
  theme(plot.title = element_text(hjust = 0.5))

# DotPlots
DotPlot(SMC.only,
        features = c("Myh11", "Acta2", "Tagln", "Cnn1", "Myocd", "Smtn", "Cald1", "Itga8")
) +
  scale_color_viridis(direction = 1, option = "D") +
  ggtitle("Expression of contractile SMC genes") +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))

DotPlot(SMC.only,
        features = c("Fn1", "Col1a1", "Lum", "Dcn", "Bgn", "S100a4", "Mmp2")
) +
  scale_color_viridis(direction = 1, option = "D") +
  ggtitle("Expression of FBMC genes") +
  theme(axis.text.x = element_text(angle = 90), plot.title = element_text(hjust = 0.5))