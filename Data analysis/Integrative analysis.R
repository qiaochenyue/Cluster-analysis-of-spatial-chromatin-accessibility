#Loom to SingleCellExperiment
https://github.com/cellgeni/sceasy
sceasy::convertFormat('filename.loom', from="loom", to="sce",
                      outFile='filename.rds')
seRNA <- readRDS("l5_all.rds")

#使用 addCoAccessibility()存储co-accessibility information
proj_in_tissue <- addCoAccessibility(
  ArchRProj = proj_in_tissue,
  reducedDims = "IterativeLSI"
)
#通过getCoAccessibility()检索信息，returnLoops = FALSE表示返回DataFrame
cA <- getCoAccessibility(
  ArchRProj = proj_in_tissue,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA #这个DataFrame中，queryHits和subjectHits列表示发现相关的两个峰值的索引。correlation列出来两个峰值之间的相关性。
metadata(cA)[[1]] #这个DataFrame中还包含一个元数据组件，其中包含相关峰值的GRanges对象。

#set returnLoops = TRUE
cA <- getCoAccessibility(
  ArchRProj = proj_in_tissue,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = TRUE
)
cA[[1]]

#Plotting browser tracks of Co-accessibility
markerGenes  <- c(
  "Nova1", 
  "Nova2",
  "Rbfox3", 
  "Srrm4",
  "Tbr1", "Satb2", "Zbtb20"
)

p <- plotBrowserTrack(
  ArchRProj = proj_in_tissue, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(proj_in_tissue)
)
grid::grid.newpage()
grid::grid.draw(p$Nova2)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
        ArchRProj = proj_in_tissue, 
        addDOC = FALSE, width = 5, height = 5)

#To identify peak-to-gene links in ArchR, we use the addPeak2GeneLinks() function.
projCUTA <- addPeak2GeneLinks(
  ArchRProj = projCUTA,
  reducedDims = "IterativeLSI"
)
#
p2g <- getPeak2GeneLinks(
  ArchRProj = projCUTA,
  corCutOff = 0.45,
  resolution = 1,
  retu
  returnLoops = FALSE
)

p2g <- getPeak2GeneLinks(
  ArchRProj = projCUTA,
  corCutOff = 0.45,
  resolution = 1, #(1,1000,10000)
  returnLoops = TRUE
)
markerGenes  <- c(
  "Nova1", 
  "Nova2",
  "Rbfox3", 
  "Srrm4",
  "Tbr1", "Satb2", "Zbtb20"
)
p <- plotBrowserTrack(
  ArchRProj = projCUTA, 
  groupBy = "predictedGroup", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(projCUTA)
)
grid::grid.newpage()
grid::grid.draw(p$CD14)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
        ArchRProj =proj_in_tissue, 
        addDOC = FALSE, width = 5, height = 5)
metadata(p2g)[[1]]