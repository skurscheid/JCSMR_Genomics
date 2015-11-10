# DANPOS2 post-processing

#------------load libraries------------------
library("GenomicFeatures")
library("ChIPpeakAnno")
library("ggplot2")
library("Gviz")
library("GenomicRanges")
library("rtracklayer")
library("GenomicFeatures")

#------------import external functions------------------
source("~/Development/GeneralPurpose/R/heatmap.3.R")

danpos2.results <- read.table("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT/Users_u1001407_Data_Tremethick_EMT_GenomeWide_Conditions_TGFb_ChIP-Users_u1001407_Data_Tremethick_EMT_GenomeWide_Conditions_WT_ChIP.positions.integrative.xls",
                              header = T, 
                              as.is = T,
                              sep = "\t")

gr.danpos2.results <- GRanges(danpos2.results$chr, IRanges(danpos2.results$start, danpos2.results$end), strand = "*", danpos2.results[, c(4:23)])
gr.TGFb_H2AZ_ChIP_bgsub_Fnor <- import("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT/pooled/TGFb_ChIP.bgsub.Fnor.smooth.bw") 
gr.WT_H2AZ_ChIP_bgsub_Fnor <- import("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT/pooled/WT_ChIP.bgsub.Fnor.smooth.bw")
gr.TGFb_vs_WT_diff <- import("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/TGFb_vs_WT/diff/TGFb_vs_WT.pois_diff.bw")

heatmap.3(mcols(subsetByOverlaps(gr.danpos2.results, gr.MSigDB.EMT_associated.cfam.tss1500))[,c("control_smt_val", "treat_smt_val")])

# create a histogram of the complete data (here log2 transformed)
df1 <- as(log2(danpos2.results[, c("control_smt_val")] + 0.0001), "matrix")
df1 <- rbind(df1, as(log2(danpos2.results[, c("treat_smt_val")] + 0.0001), "matrix"))
df1 <- data.frame(df1)
df1$var <- c(rep("ctrl", nrow(danpos2.results)), rep("treat", nrow(danpos2.results)))
hm1 <- ggplot(df1,aes(x=df1, group=var))
hm1 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))

# histogram of summit counts of nucleosomes located in TSS+/-1500 of TGFb-induced EMT genes
df2 <- as(mcols(subsetByOverlaps(gr.danpos2.results, gr.MSigDB.TGFb_induced_EMT.cfam.tss1500))[,c("control_smt_val", "treat_smt_val")], "data.frame")
rownames(df2) <- mcols(subsetByOverlaps(gr.danpos2.results, gr.MSigDB.TGFb_induced_EMT.cfam.tss1500))$row_id
heatmap.3(as.matrix(log2(df2 + 1)), trace = "none")


# Using Gviz for visualization of some of the data
i <- 2
gr1 <- promoters(gr.mesenchymalMarkers.genes[i], upstream = 10000, downstream = 10000)
dT.WT <- DataTrack(subsetByOverlaps(gr.WT_H2AZ_ChIP_bgsub_Fnor, gr1), type = "h", col = "black", name = "Control")
dT.TGFb <- DataTrack(subsetByOverlaps(gr.TGFb_H2AZ_ChIP_bgsub_Fnor, gr1), type = "h", col = "grey", name = "TGFb")
dT.Diff <- DataTrack(subsetByOverlaps(gr.TGFb_vs_WT_diff, gr1), type = "h", col = "blue", name = "Difference [+/- log10p-val]")

biomTrack <- BiomartGeneRegionTrack(genome = "canFam3", 
                                    chromosome = as(seqlevels(gr1)[i], "character"),
                                    start = as(start(gr1), "integer"),
                                    end = as(end(gr1), "integer"),
                                    name = paste(mcols(gr1)$hgnc_symbol),
                                    mart = dog)
gat <- GenomeAxisTrack()
plotTracks(list(gat, biomTrack, dT.WT, dT.TGFb, dT.Diff))


# trying to directly visualize the DANPOS2 resuls from the integrative presentation of data
DT1 <- DataTrack(subsetByOverlaps(gr.danpos2.results, gr1), data = mcols(subsetByOverlaps(gr.danpos2.results, gr1))$control_smt_val, type = "l")
DT2 <- DataTrack(subsetByOverlaps(gr.danpos2.results, gr1), data = mcols(subsetByOverlaps(gr.danpos2.results, gr1))$treat_smt_val, type = "l")
DT3 <- DataTrack(subsetByOverlaps(gr.danpos2.results, gr1), data = mcols(subsetByOverlaps(gr.danpos2.results, gr1))$smt_log2FC,type = "l")
DT4 <- DataTrack(subsetByOverlaps(gr.danpos2.results, gr1), data = -1 * log10(mcols(subsetByOverlaps(gr.danpos2.results, gr1))$smt_diff_FDR), type = c("p", "g"))

max.y <- max(max(values(DT1)), max(values(DT2)))
displayPars(DT1) <- list(ylim = c(0,max.y))
displayPars(DT2) <- list(ylim = c(0,max.y))

  
plotTracks(list(gat, biomTrack, DT1, DT2, DT3, DT4), from = start(gr1), to = end(gr1))

df1 <- as(mcols(subsetByOverlaps(gr.danpos2.results, gr1[1]))$control_smt_val, "matrix")
df1 <- rbind(df1, as(mcols(subsetByOverlaps(gr.danpos2.results, gr1[1]))$treat_smt_val, "matrix"))
df1 <- data.frame(df1)
df1$var <- c(rep("ctrl", length(subsetByOverlaps(gr.danpos2.results, gr1[1]))), 
             rep("treat", length(subsetByOverlaps(gr.danpos2.results, gr1[1]))))
df1$pos <- rep(c(1:length(subsetByOverlaps(gr.danpos2.results, gr1[1]))), 2)
p <- ggplot(df1, aes(x = pos, y = df1 , group = var, colour = var))

#------------prepare annotation data------------------
# create Canis familiaris TXDB object for peak annotation
chromInfo <- read.table("/Volumes/gduserv/Data/RefGenomes/Canis_familiaris/Ensembl/chromInfo.txt", header = F, as.is = T, sep = "\t")
colnames(chromInfo) <- c("chrom", "length")
TxDb.Cfam3.Ensembl <- makeTxDbFromGFF("/Volumes/gduserv/Data/Annotations/CanFam3/Canis_familiaris.CanFam3.1.82.chr.gtf", 
                                      organism = "Canis familiaris", 
                                      chrominfo = chromInfo)
Cfam3.genes <- genes(TxDb.Cfam3.Ensembl)

#------------annotate peaks------------------ 
danpos2.anno <- annotatePeakInBatch(gr.danpos2.results, AnnotationData = Cfam3.genes)

# only consider peaks/nucleosome position upstream or on the annotated TSS
upTSS <- which(mcols(danpos2.anno)$insideFeature %in% c("overlapStart", "upstream"))
gr.danpos2.upTSS <- danpos2.anno[upTSS]
hist(mcols(gr.danpos2.upTSS)$smt_diff_FDR)
gr.danpos2.upTSS <- gr.danpos2.upTSS[which(mcols(gr.danpos2.upTSS)$smt_diff_FDR < 0.001)]
gr.danpos2.upTSS.5kb <- gr.danpos2.upTSS[which(mcols(gr.danpos2.upTSS)$distancetoFeature <= 5000)]

df3 <- as(log2(mcols(gr.danpos2.upTSS)$"control_smt_val" + 1), "matrix")
df3 <- rbind(df3, as(log2(mcols(gr.danpos2.upTSS)$"treat_smt_val" + 1), "matrix"))
df3 <- data.frame(df3)
df3$var <- c(rep("ctrl", length(gr.danpos2.upTSS)), rep("treat", length(gr.danpos2.upTSS)))
hm2 <- ggplot(df3,aes(x=df3, group=var))
hm2 + geom_histogram(alpha = 0.6, position = "identity", group= var)
hm2 + geom_histogram(alpha = 0.6, position = "identity", aes(y = ..density..)) + geom_density(alpha = 0.4, position = "identity", aes(color = var))

df3.1 <- mcols(gr.danpos2.upTSS)[c("treat_smt_val", "control_smt_val")]
df3.1 <- as.matrix(df3.1)
heatmap.3(log2(df3.1 + 1), trace = "none", cexCol = 0.7)

# check GO enrichment of peaks
GO.danpos2.upTSS.5kbr <- getEnrichedGO(gr.danpos2.upTSS.5kb, orgAnn = "org.Cf.eg.db", maxP=0.1, multiAdj = T, minGOterm = 10, multiAdjMethod = "BH")
GO.danpos2.anno <- getEnrichedGO(danpos2.anno, orgAnn = "org.Cf.eg.db", maxP=0.1, multiAdj = T, minGOterm = 10, multiAdjMethod = "BH")

# look at control and treatment enriched separately
danpos2.results.subset <- 
GO.danpos2.anno_ctrl_enriched <- getEnrichedGO(danpos2.anno[which(mcols(danpos2.anno)$smt_log2FC < -5)], orgAnn = "org.Cf.eg.db", maxP=0.1, multiAdj = F, minGOterm = 10, multiAdjMethod = "BH")
GO.danpos2.anno_treat_enriched <- getEnrichedGO(danpos2.anno[which(mcols(danpos2.anno)$smt_log2FC > 5)], orgAnn = "org.Cf.eg.db", maxP=0.1, multiAdj = F, minGOterm = 10, multiAdjMethod = "BH")
reactome.danpos2.anno_ctrl_enriched <- getEnrichedPATH(danpos2.anno[which(mcols(danpos2.anno)$smt_log2FC < -5)], orgAnn = "org.Cf.eg.db", pathAnn = "reactome.db", maxP = 0.1, minPATHterm = 10)
reactome.danpos2.anno_treat_enriched <- getEnrichedPATH(danpos2.anno[which(mcols(danpos2.anno)$smt_log2FC > 5)], orgAnn = "org.Cf.eg.db", pathAnn = "reactome.db", maxP = 0.1, minPATHterm = 10)



