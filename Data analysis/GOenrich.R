write.csv(markerList_pos$C1,file = "ME11_C1_markers.csv")

library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)

## Read markers genes
data_dir <- './ME11_C1_markers.txt'
markerList <- read.table(data_dir, header = TRUE, stringsAsFactors = FALSE)
markers <- markerList$name


data(geneList, package = "DOSE")
g_list <- names(geneList)[1:100]
head(g_list)
[1] "4312"  "8318"  "10874" "55143" "55388" "991"  
#从DOSE包里取出了100个基因的EntrezID。


eG <- enrichGO(gene = g_list, #需要分析的基因的EntrezID
               OrgDb = org.Hs.eg.db, #人基因数据库
               pvalueCutoff =0.01, #设置pvalue界值
               qvalueCutoff = 0.01, #设置qvalue界值(FDR校正后的p值）
               ont="all", #选择功能富集的类型，可选BP、MF、CC，这里选择all。
               readable =T)
#富集分析的类型可选BP（biological process）、MF（molecular function）、CC（Cellular Component）。

#输出结果为txt文件
write.table(eG,file="eG.txt", sep="\t", quote=F, row.names = F)

#barplot
pdf(file="eGO_barplot.pdf",width = 8,height = 10) 
barplot(eG, x = "GeneRatio", color = "p.adjust", #默认参数（x和color可以根据eG里面的内容更改）
        showCategory =10, #只显示前10
        split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free') #以ONTOLOGY类型分开绘图
dev.off()
#dotplot
dotplot(ego,x = "GeneRatio", color = "p.adjust", size = "Count", 
        showCategory =5,
        split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free') 