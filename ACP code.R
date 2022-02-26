library(Seurat)
library(SingleR)
library(DoubletFinder)
library(tidyverse)
library(patchwork)
library(dplyr)
library(GSVA)
library(infercnv)
library(AnnoProbe)
library(ComplexHeatmap)
library(limma)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(monocle3)
library(SCENIC)


## Create SeuratObject
ACP1.data<-Read10X("CRA-ACP-1")
ACP2.data<-Read10X("CRA-ACP-2")
ACP3.data<-Read10X("CRA-ACP-3")
ACP1 <- CreateSeuratObject(ACP1.data, project = "ACP1", min.cells = 3, min.features = 200)
ACP2 <- CreateSeuratObject(ACP2.data, project = "ACP2", min.cells = 3, min.features = 200)
ACP3 <- CreateSeuratObject(ACP3.data, project = "ACP3", min.cells = 3, min.features = 200)

## Quality Control
# Calculate mitochondrial gene ratios
ACP1[["percent.mt"]] <- PercentageFeatureSet(ACP1, pattern = "^MT-")
# Calculate ribosomal gene ratios
ACP1[["percent.rb"]] <- PercentageFeatureSet(ACP1, pattern = "^RP[LS]")
# Calculate erythrocyte-related gene ratios
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(ACP1))
ACP1[["percent.HB"]]<-PercentageFeatureSet(ACP1, features=HB.genes) 

## Check quality control indicators
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB")
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(ACP1, group.by=group, pt.size = 0.01,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2)    
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 9, height = 8)

## Set Cell Screening Thresholds
minGene=500
maxGene=5500
maxUMI=20000
pctMT=15
pctHB=1

ACP1 <- subset(ACP1, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & 
                 nCount_RNA < maxUMI & percent.mt < pctMT & percent.HB < pctHB)

## inferCNV package
groupinfo <- data.frame(v1=colnames(ACP1),
                        v2=seurat_clusters)
geneInfor=annoGene(rownames(ACP1),"SYMBOL",'human')
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
ACP1.expr <- ACP1@assays$RNA@counts %>% as.data.frame()
ACP1.expr <- ACP1.expr[rownames(ACP1.expr)%in% geneInfor[,1],]

write.table(ACP1.expr,file ='ACP1.expr.txt',sep = '\t',quote = F)
write.table(groupinfo,file='groupFiles.txt',sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file= 'geneFile.txt',sep = '\t',quote = F,col.names = F,row.names = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix="ACP1.expr.txt",
                                    annotations_file="groupFiles.txt",
                                    delim="\t",
                                    gene_order_file="geneFile.txt",
                                    ref_group_names=c("immune.cell"))
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1, 
                              out_dir="~/project",
                              cluster_by_groups=F, 
                              denoise=TRUE,
                              HMM=TRUE,
                              num_threads=6)

## Sample integration and clustering
ACP.merge <- merge(ACP1,ACP2,ACP3)
scRNA <- ACP.merge
scRNAname <- "ACP"
tmp <- colnames(scRNA@meta.data)
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = scRNA@meta.data[,tmp])
scRNA.list <- SplitObject(scRNA, split.by = "orig.ident")
for (i in 1:length(scRNA.list)) {
  scRNA.list[[i]] <- NormalizeData(scRNA.list[[i]], verbose = FALSE)
  scRNA.list[[i]] <- FindVariableFeatures(scRNA.list[[i]], selection.method = "vst", nfeatures = 3000,
                                          verbose = FALSE)}
## Find anchors
scRNA.anchors <- FindIntegrationAnchors(object.list = scRNA.list,
                                        normalization.method = "LogNormalize",
                                        anchor.features = 3000)
## Integrate
scRNA.integrated <- IntegrateData(anchorset = scRNA.anchors,
                                  normalization.method = "LogNormalize")
## Reduction
scRNA.integrated <- ScaleData(scRNA.integrated,verbose = FALSE,
                              vars.to.regress = c("percent.mt"))
DefaultAssay(scRNA.integrated) <- "integrated"
scRNA <- RunPCA(scRNA.integrated, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA, ndims=50)
pc.num=1:35
scRNA <- scRNA %>% RunTSNE(dims=pc.num) %>% RunUMAP(dims=pc.num)
resnum <- 0.5
scRNA <- FindNeighbors(scRNA,dims = pc.num) %>% FindClusters(resolution = resnum)
##SingelR
testdata <- GetAssayData(scRNA, slot="data")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = ref, labels = ref$label.fine, clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
scRNA@meta.data$SingleR = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA$seurat_clusters == celltype$ClusterID[i]),'SingleR'] <- celltype$celltype[i]}

ACP <- scRNA
## UMAP plot
CairoPDF("ACP.UMAP",width = 8,height = 6)
DimPlot(ACP, reduction = "umap", label = F, group.by = "seurat_clusters", pt.size = 0.5)
dev.off()
## tSNE plot
CairoPDF("ACP.tSNE",width = 8,height = 6)
DimPlot(ACP, reduction = "tsne", label = F, group.by = "seurat_clusters", pt.size = 0.5)
dev.off()

sweep.res.list <- paramSweep_v3(ACP, PCs = pc.num, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = 0.09
homotypic.prop <- modelHomotypic(ACP$seurat_clusters)
nExp_poi <- round(DoubletRate*ncol(ACP)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
ACP <- doubletFinder_v3(ACP, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                         nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

## Differential gene analysis

ACP.DEG <- FindAllMarkers(ACP.DEG, logfc.threshold = 0.25, slot = "data", only.pos = T)
write.csv(ACP.DEG, "ACP.DEG.csv")

## Heatmap showing gene expression in different clusters
ACP.expr <- ACP@assays$RNA@scale.data
col.anno <- HeatmapAnnotation(
  celltype = heatmap.metadata$celltype,
  cluster = heatmap.metadata$cluster,
  col = list(cluster = col_cluster,
             celltype = col_celltype))
row.anno <- rowAnnotation(
  row.label = anno_mark(at = which(!ACP.marker$label == "NA"),labels = ACP.marker$gene[which(!ACP.marker$label == "NA")]))

myBreaks <- seq(-3, 3, length.out = 99)
myCol <- colorRampPalette(c("#0A60A7", "white", "#B42124"))(99)

CairoPDF("ACP.heatmap",width = 10,height = 10)
ComplexHeatmap::Heatmap(
  ACP.expr,
  col = colorRamp2(myBreaks, myCol),
  cluster_rows = T,
  cluster_columns = T,
  clustering_method_rows="ward.D2",
  show_column_names = F,
  show_row_names = F,
  top_annotation = col.anno,
  column_split = heatmap.metadata$cluster,
  cluster_column_slices = F,
  column_gap = unit(0, "mm"),
  right_annotation = row.anno,
  row_split = cluster,
  cluster_row_slices = F,
  row_gap = unit(0, "mm"),
  show_column_dend=F,
  show_row_dend = F)
dev.off()

## non-negative matrix factorization
nmfdat <- ACP1@assays$RNA@scale.data
nmfdat[nmfdat < 0] <- 0
print(system.time(res_10 <- nmf(nmfdat, rank = 10, method="snmf/r", seed=123456))) 
signature <- NMF::basis(res_10)

FeatureScore <- extractFeatures(res_10, 200L)
FeatureScore <- lapply(FeatureScore, function(x) rownames(res_10)[x])
FeatureScore <- do.call("rbind", FeatureScore)
rownames(FeatureScore) <- paste0("res_10", 1:5)
DT::datatable(t(FeatureScore))

topRank <- 100
programG <- list()
for (i in 1:length(ACPlist)){
  filedir <- paste0("./NNMF.snmf/", ACPlist[i], "/signature", ".txt")
  geneloading <- read.table(filedir, header = T, sep = "\t")
  geneloading$maxC <- apply(geneloading, 1, which.max) %>% paste0(ACPlist[i], "_", .)
  
  topgenelist <- rownames_to_column(geneloading, var = "gene") %>%
    pivot_longer(., cols = starts_with("A"), 
                 names_to = "program", values_to = "loading")
  
  topgenelist <- dplyr::filter(topgenelist, maxC == program) %>% 
    group_by(maxC) %>% top_n(n = topRank, wt = loading)
  topgenelist <- split(topgenelist$gene, topgenelist$maxC)
  programG <- c(programG, topgenelist)
}

ACP.Nordata <- ACP@assays$RNA@data
cells_rankings <- AUCell_buildRankings(as.matrix(ACP.Nordata), nCores = 2, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets = programG, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
programAUC <- getAUC(cells_AUC)

M <- cor(t(programAUC), method = "pearson")
NNMF.corrplot <- corrplot(M, 
                          method = "color", 
                          order = "hclust", 
                          hclust.method = "ward.D2", 
                          tl.pos = "n",
                          col = rev(brewer.pal(n = 8, name = "RdBu"))) 

## GSVA
scRNA <- ACP
scRNAname <- "ACP"
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA)
Idents(scRNA) <- "redefine_clusters"
expr <- AverageExpression(scRNA, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]
expr <- as.matrix(expr)
## msigdbr GO geneset
ref <- "GOBP"
genesets <- msigdbr(species = "Homo sapiens", category = "C5") 
genesets <- subset(genesets, gs_subcat=="GO:BP", select = c("gs_name", "gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)
## KEGG geneset
kegg.total_genesets<-read.gmt("~/Total_kegg_wangyue.gmt") %>% as.matrix() %>% as.data.frame()
kegg_name <- read_csv("~/KEGG_name.csv")
ref <- "KEGG_metabolism"
kegg.metabolism_name <- kegg_name$Metabolism
kegg.metabolism_genesets <- kegg.total_genesets[kegg.total_genesets$term %in% kegg.metabolism_name,]
kegg.metabolism_genesets <- split(kegg.metabolism_genesets$gene, kegg.metabolism_genesets$term)
genesets <- kegg.metabolism_genesets
## Hallmark geneset
ref <- "HallMark"
genesets <- msigdbr(species = "Homo sapiens", category = "H") 
genesets <- subset(genesets, select = c("gs_name","gene_symbol")) %>% as.data.frame()
genesets <- split(genesets$gene_symbol, genesets$gs_name)

## run gsva
gsvamethod <- "gsva"
gsva.res <- gsva(expr, genesets, method= gsvamethod)
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)

## Identification of differential pathways by limma
scRNA.metadata <- scRNA@meta.data %>% as.data.frame()
df.group <- data.frame(cell=rownames(scRNA.metadata),
                       group=scRNA.metadata$sampletype,
                       stringsAsFactors = F)
cluster.group <- df.group$group
cluster.group[which(cluster.group=="ACP")] <- 1
cluster.group[which(cluster.group=="control")] <- 0
cluster.group <- as.numeric(cluster.group)
design <- cbind(sampleGroup1=1, sampleGroup2vs1=cluster.group)
fit <- lmFit(gsva.res, design)
fit <- eBayes(fit)
sigPathways <- topTable(fit, coef="sampleGroup2vs1", 
                        number=Inf, p.value=1, adjust="BH")


## Monocle3
mono3.scRNA <- scRNA
data <- GetAssayData(mono3.scRNA, assay = 'RNA', slot = 'counts')
cell_metadata <- mono3.scRNA@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('cds.umap')
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(mono3.scRNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
k.num <- 20
cds <- cluster_cells(cds,k = k.num)
cds <- learn_graph(cds,close_loop=T)
plot_cells(cds,
           label_groups_by_cluster = F,
           label_leaves = F,
           label_branch_points = F,
           label_roots = F,
           color_cells_by = "celltype",
           label_cell_groups = F,
           cell_size = 0.7)
rownames(cds@principal_graph_aux[["UMAP"]]$dp_mst) <- NULL
colnames(cds@int_colData@listData$reducedDims@listData$UMAP) <- NULL
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime",
           label_cell_groups = F,
           label_leaves = F,
           label_branch_points = F,
           label_roots = F,
           show_trajectory_graph = F,
           cell_size = 0.7)
Track_genes_sig <- c("LEF1","MUC1")
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="mono3.celltype",
                         min_expr = min_exprlevel)

## Addmodulescore
Signature_genelist <- geneset
ACP <- AddModuleScore(ACP,
                      features = list(Signature_genelist),
                      ctrl = 50,
                      nbin = 24,
                      name = "xxScore")

## SCENIC analysis
exprMat <- GetAssayData(scRNA, assay = 'RNA', slot = 'data') %>% as.matrix()
cellInfo <- data.frame(scRNA@meta.data) %>% dplyr::select(celltype, seurat_clusters)
mydbDIR <- "~/project/SCENIC_ref_hg38/"
mydbs <- dir(mydbDIR)
names(mydbs) <- c("10kb", "500bp")

scenicOptions <- initializeScenic(org = "hgnc", 
                                  nCores = 6,
                                  dbDir = mydbDIR, 
                                  dbs = mydbs,
                                  datasetTitle = "HNC")
genesKept <- geneFiltering(exprMat, scenicOptions, 
                           minCountsPerGene = 0.015 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.015)
exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

runCorrelation(exprMat_filtered, scenicOptions)
runGenie3(exprMat_filtered, scenicOptions, nParts = 20)
runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions@settings$nCores <- 2
runSCENIC_2_createRegulons(scenicOptions)
scenicOptions@settings$nCores <- 2
runSCENIC_3_scoreCells(scenicOptions, exprMat=exprMat)
runSCENIC_4_aucell_binarize(scenicOptions, exprMat=exprMat)

Regulon <- readRDS("int/3.4_regulonAUC.Rds")
Regulon <- Regulon@assays@data@listData$AUC
Regulon_all <- Regulon[,colnames(scRNA)]
Regulon_high <- Regulon[!grepl("extended", rownames(Regulon)),]
Regulon_low <- Regulon[grepl("extended", rownames(Regulon)),]
scRNA[["Regulon"]] <- CreateAssayObject(counts = Regulon_all)
scRNA@assays$Regulon@data <- scRNA@assays$Regulon@counts

binRegulon <- readRDS("int/4.1_binaryRegulonActivity.Rds")
binRegulon_all <- binRegulon[,colnames(scRNA)]
binRegulon_high <- binRegulon[!grepl("extended", rownames(binRegulon)),]
binRegulon_low <- binRegulon[grepl("extended", rownames(binRegulon)),]
scRNA[["binRegulon"]] <- CreateAssayObject(counts = binRegulon_all)
scRNA@assays$binRegulon@data <- scRNA@assays$binRegulon@counts











