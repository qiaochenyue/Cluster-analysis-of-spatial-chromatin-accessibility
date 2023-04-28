#显示cluster
p1 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
#显示细胞类型
p2 <- plotEmbedding(ArchRProj = projHeme5, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
#出图
ggAlignPlots(p1,p2,type = "h")

#创建轨迹的第一步是以细胞组标签的有序向量的形式创建轨迹主干
trajectory <- c("Radial glia", "Postmitotic premature neurons", "Excitatory neurons")
trajectory2 <- c("Radial glia", "Oligodendrocyte Progenitors", "Premature oligodendrocyte")
trajectory3 <- c("Radial glia", "Inhibitory neuron progenitors", "Inhibitory interneurons")

trajectory <- c("Neural progenitor cells","Radial glia")

#使用addTrajectory()函数创建一条轨迹，并将其添加到ArchRProject中。使用name为轨迹命名，这样做可以在cellColData中创建一个名为"Neuron_U"的新列，它储存了轨迹中每个细胞的拟时间值。不属于轨迹的细胞用NA标记。
projCUTA <- addTrajectory(
  ArchRProj = proj_in_tissue, 
  name = "Neuron_U", 
  groupBy = "predictedGroup",
  trajectory = trajectory, 
  embedding = "UMAP", 
  force = TRUE
)
#我们可以查看这些信息，发现每个细胞都有一个在0-100之间的唯一伪时间值，并排除NA值的细胞。
head(projCUTA$Neuron_U[!is.na(projCUTA$Neuron_U)])

#使用plotTrajectory()函数绘制轨迹，该函数覆盖了UMAP嵌入的拟时间值，并显示了一个箭头，接近拟合后的轨迹路径。不属于轨迹的细胞呈灰色。trajectory告诉ArchR哪个细胞子集我们感兴趣。使用colorBy = "cellColData"告诉ArchR在cellColData中查找由name指定的列--"Neuron_U"，以此给细胞子集上色。
p <- plotTrajectory(ArchRProj = projCUTA,trajectory = "Neuron_U", colorBy = "cellColData", name = "Neuron_U")
p[[1]]
plotPDF(p,name = "Neuron_U.pdf",ArchRProj = projCUTA,addDOC = FALSE,width = 5,height = 5)

#将impute weights添加到projCUTA中
projCUTA <- addImputeWeights(projCUTA)
#绘制出"Neuron_U"轨迹，但通过Nova2基因的基因评分值给细胞着色。colorBy表示要使用的矩阵，name表示要使用的特征。
p1 <- plotTrajectory(projCUTA,trajectory = "Neuron_U", colorBy = "GeneScoreMatrix",name = "Nova2",continuousSet = "horizonExtra")
p2 <- plotTrajectory(projCUTA,trajectory = "Neuron_U", colorBy = "GeneIntegrationMatrix",name = "Nova2",continuousSet = "blueYellow")
ggAlignPlots(p1[[1]],p2[[1]],type ="h")

trajMM <- getTrajectory(ArchRProj = projCUTA, name = "Neuron_U", useMatrix = "MotifMatrix", log2Norm = FALSE)
#将SummarizedExperiment 传递给plotTrajectoryHeatmap()函数
p1 <- plotTrajectoryHeatmap(trajMM,pal = paletteContinuous(set = "solarExtra"))

trajGSM <- getTrajectory(ArchRProj = projCUTA, name = "Neuron_U", useMatrix = "GeneScoreMatrix", log2Norm = FALSE)
p2 <- plotTrajectoryHeatmap(trajGSM,pal = paletteContinuous(set = "horizonExtra"))

trajGIM <- getTrajectory(ArchRProj = projCUTA, name = "Neuron_U", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,pal = paletteContinuous(set = "blueYellow"))

trajPM <- getTrajectory(ArchRProj = projCUTA, name = "Neuron_U", useMatrix = "PeakMatrix", log2Norm = FALSE)
p4 <- plotTrajectoryHeatmap(trajPM,pal = paletteContinuous(set = "solarExtra"))

plotPDF(p1,p2,p3,p4,name = "Neuron_U-Heatmaps.pdf",ArchRProj = projCUTA, addDOC = FALSE, width = 6,height = 8)

#先找出伪时间分布中与TF基因评分相关的motif可接近性。 correlateTrajectories()输出的是一个列表对象，其中包含一个DataFrame 对象作为列表中的第一个条目。这个DataFrame包含名为idx1,matchname1,name1和VarAssay1的列，它们对应于传递给correlateTrajectories()函数的第一个轨迹（基因得分）中的索引，匹配名称，为修改名称和特征方差分位。“方差分位数”是给定特征的标准化度量，它允许我们在不同的分析之间获得相关性。该DataFrame包含了满足correlateTrajectories()函数中指定的cutoffs的所有特性。
corGSM_MM <- correlateTrajectories(trajGSM,trajMM)
#将对应的轨迹SummarizedExperiment对象进行子集化，使其只包含上面传递的显著性元素。
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1,]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2,]
#为了最好排列这些特征，我们可以创建一个新的轨迹，其中这两个轨迹的值相乘。这将允许我们创建并拍的热图，并按照“行”进行相同的排序。
trajCombined <- trajGSM2
assay(trajCombined,withDimnames=FALSE) <- t(apply(assay(trajGSM2),1,scale)) + t(apply(assay(trajMM,2),1,scale))
#从plotTrajectoryHeatmap()函数的返回值中提取最优行序
combinedMat <- plotTrajectoryHeatmap(trajCombined,returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat),rownames(trajGSM2))
#创建成对的热图
#创建基因评分轨迹的热图。通过rowOrdwe参数指定所需的行顺序
ht1 <- plotTrajectoryHeatmap(trajGSM2,pal = paletteContinuous(set = "horizonExtra"),varCutoff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2,pal = paletteContinuous(set = "solarExtra"),varCutoff = 0, rowOrder = rowOrder)

#使用SpatialPlot_traj.R包
meta.data.integration <- as.data.frame(getCellColData(ArchRProj = proj_in_tissue))[, c('Neuron_U'), drop=FALSE]
new_row_names <- row.names(meta.data.integration)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data.integration) <- new_row_names
all(row.names(meta.data.integration) == colnames(spatial.obj))

spatial.obj <- AddMetaData(object = spatial.obj, metadata = meta.data.integration)

p <- SpatialPlot_traj(spatial.obj, features = "Neuron_U",  pt.size.factor = 4, image.alpha = 0, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22)
p