suppressMessages(library(ggplot2))
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

set.seed(42)
addArchRThreads(threads = 64)

setwd("/home/hanxue/lab/scATAC/BCY_ATAC_R/joint")
proj <- loadArchRProject("./ArchRProject", showLogo = FALSE)

# ==================== 1. Motif 富集分析 ====================
proj <- addMotifAnnotations(
    ArchRProj = proj, 
    motifSet = "cisbp", 
    name = "Motif"
    )

markersPK <- readRDS("./ArchRProject/markersPK_majorType.rds")

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPK,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

enrichMotifs

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
p1 <- ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(p1, name = "cisbp_Motif_Enrichment.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 15)

# ==================== 2. EncodeTFBS(转录因子结合位点) 富集分析 ====================
proj <- addArchRAnnotations(
    ArchRProj = proj, 
    collection = "EncodeTFBS"
)

enrichEncode <- peakAnnoEnrichment(
    seMarker = markersPK,
    ArchRProj = proj,
    peakAnnotation = "EncodeTFBS",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

enrichEncode

heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
p2 <- ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(p2, name = "EncodeTFBS_Enrichment.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 15)

# ==================== 3. bulk ATAC 富集分析 ====================
proj <- addArchRAnnotations(ArchRProj = proj, collection = "ATAC")

enrichATAC <- peakAnnoEnrichment(
    seMarker = markersPK,
    ArchRProj = proj,
    peakAnnotation = "ATAC",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

enrichATAC

heatmapATAC <- plotEnrichHeatmap(enrichATAC, n = 7, transpose = TRUE)
p3 <- ComplexHeatmap::draw(heatmapATAC, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(p3, name = "ATAC_Enrichment.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 15)

# ==================== 4. Codex TFBS 富集分析 ====================
proj <- addArchRAnnotations(ArchRProj = proj, collection = "Codex")

enrichCodex <- peakAnnoEnrichment(
    seMarker = markersPK,
    ArchRProj = proj,
    peakAnnotation = "Codex",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

enrichCodex

heatmapCodex <- plotEnrichHeatmap(enrichCodex, n = 7, transpose = TRUE)
p4 <- ComplexHeatmap::draw(heatmapCodex, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(p4, name = "Codex_Enrichment.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 15)

## 保存结果
saveArchRProject(proj, load = FALSE)
