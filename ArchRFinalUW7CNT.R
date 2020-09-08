library(ArchR)
library(Cairo)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2018)
library(TFBSTools)
library(patchwork)
library(chromVAR)
library(pheatmap)
library(monocle3)
set.seed(1)

SFtheme<-theme_bw(base_size=14) + 
  theme(panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"), 
        legend.key = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))
theme_set(SFtheme)

setwd('/mnt/c/Users/Anoop/Fred Hutchinson Cancer Research Center/scCut&Tag - General/Manuscript/code_repo/gbm/app_analysis/')
fragment.path<-'/mnt/c/Users/Anoop/Fred Hutchinson Cancer Research Center/scCut&Tag - General/Manuscript/code_repo/gbm/UW7_raw_data/fragments.tsv.gz'
sample.name='UW7.H3K27me3'
addArchRThreads(threads = 8)
addArchRGenome("hg38")
inputFiles<-fragment.path
names(inputFiles)<-'UW7.H3K27me3'

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 0, 
  filterFrags = 1000, 
  minFrags=400,
  maxFrags=20000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  TileMatParams=list(tileSize=5000),
  force=T
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP",
  LSIMethod = 1
) 

proj.prelim <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = 'CNT7_Prelim',
  copyArrows = TRUE
)

proj.prelim <- filterDoublets(ArchRProj = proj.prelim)
proj.prelim <- addIterativeLSI(ArchRProj = proj.prelim, useMatrix = "TileMatrix", iterations=1,dimsToUse = 1:50, name = "IterativeLSI",force=TRUE)
proj.prelim <- addClusters(input = proj.prelim, reducedDims = "IterativeLSI",res=1.0,force=T)
proj.prelim <- addUMAP(ArchRProj = proj.prelim, reducedDims = "IterativeLSI",force=T)
umap.prelim <- plotEmbedding(ArchRProj = proj.prelim, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(umap.prelim,name='umap.prelim.pdf')
proj.prelim<-saveArchRProject(proj.prelim)
#markers based on gene activity scores
markersGS <- getMarkerFeatures(proj.prelim, useMatrix = "GeneScoreMatrix", groupBy = "Clusters", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerListGS <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.5")
markerGenes  <- c("PTPRZ1","PTPRC","MOBP","GFAP","RBFOX3")

heatmapGS.prelim <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 2", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

heatmapGSneg.prelim <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC <= -2", 
  labelMarkers = markerGenes,
  transpose = TRUE,
  invert = TRUE
)

heatmapGS.prelim
plotPDF(heatmapGS.prelim,name='heatmapGS.prelim.pdf')

heatmapGSneg.prelim
plotPDF(heatmapGSneg.prelim,name='heatmapGS.neg.prelim.pdf')

#no genes enriched in cluster 5... decide to remove it. 
#cells.clusters<-as.data.frame(getCellColData(proj.prelim,select='Clusters'))
#keep.cells<-rownames(dplyr::filter(cells.clusters, Clusters %ni% c('C5')))
#proj<-subsetArchRProject(proj.prelim,cells=keep.cells,outputDirectory ='CNT7_Final')
proj<-proj.prelim
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", iterations=1,dimsToUse = 1:50,name = "IterativeLSI",force=TRUE)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI",force=T,res=1.0)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI",force=T)
umap <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(umap,name='umap.final.pdf')
umap.nFrags <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")
umap.doublet <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAP")
plotPDF(umap.nFrags,name='umap.final.nFrags.pdf')
proj<-saveArchRProject(proj)

cds<-exportMonocle3(
  ArchRProj = proj,
  useMatrix = 'TileMatrix',
  threads = getArchRThreads(),
  verbose = TRUE,
  binarize = T,
  logFile = createLogFile("exportMonocle3")
)

cds <- cluster_cells(cds,method='UMAP',k=20,partition_qval = 0.05,res=0.02,random_seed=1)
cds <- learn_graph(cds)
m3.traj.plot<-plot_cells(cds, color_cells_by = "cluster",show_trajectory_graph = T)
plotPDF(m3.traj.plot,name='UW7.H3K27me3.m3traj.pdf')
plot_cells(cds, color_cells_by = "cluster",show_trajectory_graph=F)
cds <- order_cells(cds)
m3.pseudotime<-plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = T)
plotPDF(m3.pseudotime,name='UW7.H3K27me3.pseudotime.pdf')
pData(cds)$pseudotime<-pseudotime(cds)
pt.matrix<-(pseudotime(cds))
colData(cds)$celltype <- as.character(clusters(cds))
colData(cds)$celltype = dplyr::recode(colData(cds)$celltype,
                                                "1"="Tumor1",
                                                "2"="Tumor2",
                                                "3"="Monocytes",
                                                "4"="Neurons",
                                                "5"="Tumor3",
                                                "6"="TumorStem",
                                                "7"="Monocytes",
                                                "8"="Oligodendrocytes",
                                                "9"="Tumor4")
clust.matrix<-as.character(clusters(cds))
celltype.matrix<-as.character(colData(cds)$celltype)
cells<-colnames(cds)
proj<-addCellColData(proj,cells=cells,data=pt.matrix,name = 'pseudotime',force=T)
proj<-addCellColData(proj,cells=cells,data=clust.matrix,name = 'm3_clust',force=T)
proj<-addCellColData(proj,cells=cells,data=celltype.matrix,name = 'celltype',force=T)
m3_clust_umap <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "m3_clust", embedding = "UMAP")
umap.pseudotime <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "pseudotime", embedding = "UMAP")
m3_clust_celltype_umap <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "celltype", embedding = "UMAP")
plotPDF(umap.pseudotime,name='UW7.H3K27me3.umap.pseudotime.pdf')
plotPDF(m3_clust_umap,name='UW7.H3K27me3.umap.m3_clust.pdf')
plotPDF(m3_clust_celltype_umap,name='UW7.H3K27me3.umap.m3.celltype.pdf')

#celltype annotation by chromatin silencing score over marker genes
proj <- addImputeWeights(proj,k=4)

gene.plot <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)

plot.PTPRZ1<-gene.plot$PTPRZ1
plotPDF(plot.PTPRZ1,name='CSS.PTPRZ1.pdf')
plot.PTPRC<-gene.plot$PTPRC
plotPDF(plot.PTPRC,name='CSS.PTPRC.pdf')
plot.MOBP<-gene.plot$MOBP
plotPDF(plot.MOBP,name='CSS.MOBP.pdf')
plot.RBFOX3<-gene.plot$RBFOX3
plotPDF(plot.RBFOX3,name='CSS.RBFOX3.pdf')


#celltype annotation using bulk chip-seq projection
load('/mnt/c/Users/Anoop/Fred Hutchinson Cancer Research Center/scCut&Tag - General/Manuscript/code_repo/gbm/Bulk_data/H3K27me3.celltypes.RData')
celltypes<-celltype.data[,c(1,5)]
projection.celltypes<-projectBulkATAC(
  ArchRProj = proj,
  seATAC = celltypes,
  reducedDims = "IterativeLSI",
  embedding = "UMAP",
  n = 250,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("projectBulkATAC")
)

p<-as.data.frame(do.call(rbind, projection.celltypes[2:1]))
table(p$Type)
p$Type[grepl("scATAC",p$Type)]<-'UW7.H3K27me3'
table(p$Type)
celltype.projection<-ggplot(p, aes(x=UMAP1, y=UMAP2, color=Type))+geom_point(size=0.5, alpha=0.4)
celltype.projection+theme()
c.p<-celltype.projection+theme()
plotPDF(c.p,name='celltype.projection.monocyte.uw7gsc.pdf')

#Project bulk neural stem cell data
load('/mnt/c/Users/Anoop/Fred Hutchinson Cancer Research Center/scCut&Tag - General/Manuscript/code_repo/gbm/Bulk_data/RSE.celltypes.nsc.RData')
projection.nsc<-projectBulkATAC(
  ArchRProj = proj,
  seATAC = nsc,
  reducedDims = "IterativeLSI",
  embedding = "UMAP",
  n = 250,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("projectBulkATAC")
)

p.nsc<-as.data.frame(do.call(rbind, projection.nsc[2:1]))
table(p$Type)
p.nsc$Type[grepl("scATAC",p.nsc$Type)]<-'UW7.H3K27me3'
table(p.nsc$Type)
celltype.projection<-ggplot(p.nsc, aes(x=UMAP1, y=UMAP2, color=Type))+geom_point(size=0.5, alpha=0.4)
celltype.projection+theme()
c.p.nsc<-celltype.projection+theme()
plotPDF(c.p.nsc,name='celltype.projection.nsc.pdf')

#Create recurrent tumor object
inputFiles.recurrent<-'/mnt/c/Users/Anoop/Fred Hutchinson Cancer Research Center/scCut&Tag - General/Manuscript/code_repo/gbm/UW7R_raw_data/fragments.tsv.gz'
names(inputFiles.recurrent)<-'UW7R.H3K27me3'
ArrowFiles.recurrent <- createArrowFiles(
  inputFiles = inputFiles.recurrent,
  sampleNames = names(inputFiles.recurrent),
  filterTSS = 0, 
  filterFrags = 1000, 
  minFrags=400,
  maxFrags=100000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  TileMatParams=list(tileSize=5000),
  force=T
)

doubScores <- addDoubletScores(
  input = ArrowFiles.recurrent,
  k = 10, 
  knnMethod = "UMAP",
  LSIMethod = 1
) 

projR <- ArchRProject(
  ArrowFiles = ArrowFiles.recurrent, 
  outputDirectory = 'CNT7R_Final',
  copyArrows = TRUE
)

projR <- filterDoublets(ArchRProj = projR)
projR <- addIterativeLSI(ArchRProj = projR, useMatrix = "TileMatrix", iterations=1,name = "IterativeLSI",force=TRUE)
projR <- addClusters(input = projR, reducedDims = "IterativeLSI",force=T)
projR <- addUMAP(ArchRProj = projR, reducedDims = "IterativeLSI",force=T)
projR <- addImputeWeights(projR,k=20)
umap.R <- plotEmbedding(ArchRProj = projR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(umap.R,name='UW7R.H3K27me3.umap.pdf')
tmR<-getMatrixFromArrow(ArrowFiles.recurrent, useMatrix = "TileMatrix", binarize = T)

gene.plot.r <- plotEmbedding(
  ArchRProj = projR, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projR)
)

plot.PTPRZ1.r<-gene.plot.r$PTPRZ1
plotPDF(plot.PTPRZ1.r,name='UW7R.CSS.PTPRZ1.pdf')

tmR@assays@data$counts<-tmR@assays@data$TileMatrix
tmR@assays@data$TileMatrix<-NULL
window_size<-5000
rs<-ArchR:::.getRowSums(ArrowFiles = ArrowFiles, useMatrix = "TileMatrix")
rowRanges(tmR)<-GRanges(seqnames = rs$seqnames, ranges = IRanges(start = (rs$idx*window_size)-4999, end = (rs$idx*window_size)))


p.recurrent<-projectData(
  projector = proj,
  projectee = tmR,
  reducedDims = "IterativeLSI",
  embedding = "UMAP",
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("projectBulkATAC")
)

p.r<-as.data.frame(do.call(rbind, p.recurrent[2:1]))
table(p.r$Type)
p.r$Type[grepl("projector",p.r$Type)]<-'UW7.H3K27me3'
p.r$Type[grepl("CNT7R",p.r$Type)]<-'UW7R.H3K27me3'
table(p.r$Type)
recurrent.projection<-ggplot(p.r, aes(x=UMAP1, y=UMAP2, color=Type))+geom_point(size=0.5, alpha=0.4)
r.p<-recurrent.projection+theme()
plotPDF(r.p,name='recurrent.projection.pdf')

projR<-saveArchRProject(projR)
proj<-saveArchRProject(proj)