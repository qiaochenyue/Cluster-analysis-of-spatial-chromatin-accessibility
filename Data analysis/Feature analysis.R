## Create ArchRProject
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  filterTSS = 0,
  filterFrags = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000)
)
ArrowFiles

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE
)
proj

## Select pixels in tissue
meta.data <- as.data.frame(getCellColData(ArchRProj = proj))
meta.data['cellID_archr'] <- row.names(meta.data)
data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"
image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), 
                       filter.matrix = filter.matrix)
meta.data.spatial <- meta.data[row.names(image@coordinates), ]
proj_in_tissue <- proj[meta.data.spatial$cellID_archr, ]
proj_in_tissue

## Data normalization and dimensionality reduction 
proj_in_tissue <- addIterativeLSI(
  ArchRProj = proj_in_tissue,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

proj_in_tissue <- addClusters(
  input = proj_in_tissue,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  force = TRUE
)

proj_in_tissue <- addUMAP(
  ArchRProj = proj_in_tissue, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

plotEmbedding(ArchRProj = proj_in_tissue, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 1.5)

proj_in_tissue <- addImputeWeights(proj_in_tissue)

## Identify the marker genes for each cluster 
markersGS <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  testMethod = "wilcoxon"
)

markerList_pos <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")

markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}
markerGenes <- unlist(markerGenes)

#可视化标记特征
mG <- c("Dusp6","CD34","FANCA","FANCD2","Slc4a1","Nova2","Rarg","GATA1","IRF8","EBF1","MS4A1","Gata2","Ascl2","Notch1","Ntng1","Car10","Sox2","Pax6","Pou3f2","Olig2")
heatmapGS <- markerHeatmap(
  seMarker =markersGS,
  cutOff = "FDR <= 0.05 & Log2FC >= 0.25",
  labelMarkers = mG,
  transpose = TRUE,
  plotLog2FC = TRUE)
ComplexHeatmap::draw(heatmapGS)

plotPDF(heatmapGS,name = "Tonsil_heatmapGS",ArchRProj = proj_in_tissue,addDOC = FALSE)

#可视化Marker的基因
p <- plotEmbedding(
  ArchRProj = proj_in_tissue,
  colorBy = "GeneScoreMatrix",
  name = mG,
  embedding = "UMAP",
  quantCut = c(0.01,0.95),
  imputeWeights = getImputeWeights(proj_in_tissue)
)
p$SPRED2  

p2 <- lapply(p,function(x){
  x+guides(color = FALSE,fill = FALSE) +
    theme_ArchR(baseSize=6.5)+
    theme(plot.margin = unit(c(0,0,0,0),"cm"))+
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
    )
})

do.call(cowplot::plot_grid,c(list(ncol = 3),p2))

qMarkerGenes <- c("CD34","CD14","DTX1","SPRED2","CD8A")
p <- plotBrowserTrack(
  ArchRProj = proj_in_tissue,
  groupBy = "Clusters",
  geneSymbol = mG,
  upstream = 50000,
  downstream = 50000)

grid::grid.draw(p$CD14)

ArchRBrowser(proj_in_tissue)

## Call peaks
proj_in_tissue <- addGroupCoverages(ArchRProj = proj_in_tissue, groupBy = "Clusters")

pathToMacs2 <- findMacs2()

proj_in_tissue <- addReproduciblePeakSet(
  ArchRProj = proj_in_tissue, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

proj_in_tissue <- addPeakMatrix(proj_in_tissue)

getAvailableMatrices(proj_in_tissue)
getPeakSet(proj_in_tissue)

if("Motif" %ni% names(proj_in_tissue@peakAnnotation)){
  proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif", force = TRUE)
}

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj_in_tissue,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 0.1")

if("Motif" %ni% names(proj_in_tissue@peakAnnotation)){
  proj_in_tissue <- addMotifAnnotations(ArchRProj = proj_in_tissue, motifSet = "cisbp", name = "Motif", force = TRUE)
}

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.05 & Log2FC >= 0.25",
  transpose = TRUE
)

## ChromVAR Deviatons Enrichment
proj_in_tissue <- addBgdPeaks(proj_in_tissue, force = TRUE)

proj_in_tissue <- addDeviationsMatrix(
  ArchRProj = proj_in_tissue, 
  peakAnnotation = "Motif",
  force = TRUE
)

markersMotifs <- getMarkerFeatures(
  ArchRProj = proj_in_tissue, 
  useMatrix = "MotifMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames = 'z'
)

p <- plotGroups(ArchRProj = proj_in_tissue,groupBy = "predictedGroup",colorBy = "MotifMatrix",name = markersMotifs13,imputeWeights = getImputeWeights(proj_in_tissue))

p2 <- lapply(seq_along(p),function(x){
  if(x!=1){
    p[[x]] + guides(color = FALSE,fill = FALSE) +
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm")) +
      theme(
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill =FALSE) +
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})

#Footprinting分析

#Motif Footprinting
motifPositions <- getPositions(proj_in_tissue)

motifs13 <- c("Nova2","Nova1","Rbfox3","Rarg","Sox2","Mir9-3","Olig2","Notch1","Mir124a-3","Mir9-2","Map1b","Ascl1","GATA1","CEBPA","EBF1","BRN2","Sp")
markersMotifs <- unlist(lapply(motifs13,function(x) grep(x,names(motifPositions),value = TRUE)))
markersMotifs

seFoot <- getFootprints(
  ArchRProj = proj_in_tissue,
  positions = motifPositions[markersMotifs],
  groupBy = "Clusters"
)

#Footprints 对于 Tn5偏差的标准化
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj_in_tissue,
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
