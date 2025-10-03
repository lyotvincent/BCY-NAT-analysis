suppressMessages(library(ggplot2))
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

set.seed(42)
addArchRThreads(threads = 64)

setwd("/home/hanxue/lab/scATAC/BCY_ATAC_R/joint")

## 因为heatmap绘制和jupyter notebook不兼容，所以在R中运行（https://github.com/GreenleafLab/ArchR/issues/1601）

proj <- loadArchRProject("./ArchRProject", showLogo = FALSE)
PeakCellType <- readRDS("./ArchRProject/markersPK_majorType.rds")

## DAR motif enrichment 差异峰值motif富集分析
enrichMotifs <- peakAnnoEnrichment(
  seMarker = PeakCellType,
  ArchRProj = proj,
  peakAnnotation = "homer",
  cutOff = "FDR < 0.01 & Log2FC >= 1"
)
enrichMotifs
heatmapMotif <- plotEnrichHeatmap(enrichMotifs, n = 7, clusterCols = F, transpose = TRUE)
p1 <- ComplexHeatmap::draw(heatmapMotif, 
                     heatmap_legend_side = "bot", 
                     annotation_legend_side = "bot")
plotPDF(p1, name = "Motif_Enrichment.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 15)
