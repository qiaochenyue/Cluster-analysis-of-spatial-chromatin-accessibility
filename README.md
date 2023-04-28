## 染色体可及性聚类分析



染色体可及性是指真核生物基因组中DNA的物理可访问性和调控因子可结合性等重要特性。研究染色体可及性可以帮助我们深入了解基因调控机制，识别组织或细胞类型中的候选调控基因组区域，以及发现潜在的治疗靶点。目前已经发展出多种生化方法来描述染色质的可及性，但对于空间分辨染色质可及性的方法和数据却相对缺乏。最新的染色质可及性研究方法——空间ATAC-seq，可以实现在组织原位对染色质可及性进行研究。



## 数据预处理

1) Raw Fastq数据包括Raw Read 1（包含空间条形码A和条形码B）和Raw Read 2（包含基因组序列）。使用两个恒定的接头序列过滤Raw Read 1，得到New Read 1（包含基因组序列）和New Read 2（包含空间条形码A和条形码B）。将Raw Read 2改为New Read 3。

2) 将New Read 1、New Read 2和New Read 3 转化为Cell Ranger ATAC格式（10x Genomics），将新生成的fastq文件与mm10参考基因组对齐，过滤以删除重复项，并使用Cell Ranger ATAC v1.2进行计数。最后生成Fragments文件用于下游分析。Fragments 文件包含关于基因组和组织位置的片段信息（条形码A x条形码B）。

## 数据分析和可视化 

1) 将文件读入ArchR中。使用IterativeLSI 方法进行降维 (iterations = 2, resolution = 0.2, varFeatures = 25000, dimsToUse = 1:30, n.start = 10)，Harmony进行批次效应矫正。使用Seurat的FindClusters功能进行聚类，在降维空间中使用UMAP方法实现可视化。

2) 基因评分和Marker基因。

3) 使用ArchR中的addReproductiblePeakSet函数用macs2 进行call peak。使用addMotifAnnotation()函数实现差异峰的motif富集（cutOff = "FDR <= 0.05 & Log2FC >= 0.1"）。chromVAR先根据所有细胞或者样本的平均情况来计算期望开放性，然后用它来计算每个注释，每个细胞或样本的偏差，最后对开放进行纠正。之后使用fetFootprints()函数进行Footprints分析。

4) 用clusterProfiler包GO富集分析。

5) 使用FindTransferAnchors函数跨平台将spatial-ATAC-seq和scRNA-seq联系。通过整合识别的scRNA-seq细胞类型为spatial-ATAC-seq数据中的细胞着色，使用scRNA-seq信息给spatial-ATAC-seq clusters做标记。并使用SpatialDimPlot_new.R包实现空间映射。

6) 使用addTrajectory函数创建轨迹，使用plotTtrajectory函数使用默认值绘制选定基因得分沿伪时间的分布。
