#qPCR data analysis
library(HTqPCR)

# set working directory
setwd("~/Data/Tremethick/EMT/GenomeWide/danpos_analysis/")

qPCRGeneList <- readLines("../../MDCK qPCR data/genelist.txt")
qPCRGenesTab <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = qPCRGeneList, dog)

qPCRGenesTab
qPCRGeneList[which(!qPCRGeneList %in% qPCRGenesTab$hgnc_symbol)]
qPCRGeneList.missing <- c("KRT7" = "ENSCAFG00000007307",
 "LOC488207" = NULL,
 "MMP3" = "ENSCAFG00000015063",
 "OCLN" = "ENSCAFG00000007805",
 "SIP1" = "ENSCAFG00000013859" ,
"STEAP1" = "ENSCAFG00000001866",
"TCF4" = "ENSCAFG00000000140",
"TGFB1" = "ENSCAFG00000005014",
"TMEFF1" = "ENSCAFG00000002577",
"TWIST1" = "ENSCAFG00000012469", #using TWIST2 - TWIST does not seem to exist in dog genome
"LOC478215/H2AZ" = "ENSCAFG00000010615",
"HPRT1" = "ENSCAFG00000018870",
"LDHAL6B" = "ENSCAFG00000009211",
"GAPDH" = "ENSCAFG00000015077"
)
t1<- cbind(rownames(as.data.frame(qPCRGeneList.missing)), qPCRGeneList.missing)
colnames(t1) <- c("hgnc_symbol", "ensembl_gene_id")
t2 <- rbind(qPCRGenesTab, t1[,c("ensembl_gene_id", "hgnc_symbol")])
rownames(t2) <- t2$hgnc_symbol
qPCRGenesTab <- t2
rm(list = c("t1", "t2"))
qPCRGenes.control <- c("B2M", "GUSB", "HPRT1", "GAPDH", "LDHAL6B")
qPCRGenesTab[qPCRGenes.control,]$type <- "CONTROL"
qPCRGenesTab$type <- "GOI"

qPCRData.H2AZ_KD <- read.csv("../../MDCK qPCR data/Shuyi EMT PCR Array Table 2.csv", header = T, row.names = 1, as.is = T)

lapply(seq_along(1:ncol(qPCRData.H2AZ_KD)), function(x) {
Ct <- qPCRData.H2AZ_KD[,x]
genes <- rownames(qPCRData.H2AZ_KD)
pos <- c(1:nrow(qPCRData.H2AZ_KD))
flag <- rep("OK", nrow(qPCRData.H2AZ_KD))
type <- c(rep("GOI", 84), rep("HKG", 5))
fn <- paste("../../MDCK qPCR data/", colnames(qPCRData.H2AZ_KD)[x], ".csv", sep = "")
print(fn)
tab <- data.frame(cbind(genes, Ct, pos, flag, type ))
write.table(tab, file = fn, row.names = F, col.names = F, sep = "\t")
})

column.info <- list(flag = "flag", Ct = "Ct", feature = "genes", position = "pos", type = "type")

qPCRdata.raw <- readCtData(files$File, 
                           format = "plain", 
                           path = "../../MDCK qPCR data", 
                           n.features = 89, 
                           column.info = list(flag = 4, Ct = 2, feature = 1, position = 3, type = 5), 
                           n.data = 1,
                           header = F)

pData(qPCRdata.raw) <- data.frame(Sample = c("Control", "Control", "H2AZ_KD", "H2AZ_KD"), Replicate = rep(1:2,2 ))


g <- featureNames(rawData)[1:10]
plotCtOverview(qPCRdata.raw, calibrator = "Control", genes = g, conf.int = T, ylim = c(0, 2), groups = c("Control", "Control", "H2AZ_KD", "H2AZ_KD"))

sr.norm <- normalizeCtData(qPCRdata.raw, norm = "scale.rank")
q.norm <- normalizeCtData(qPCRdata.raw, norm = "quantile")
dCT.norm <- normalizeCtData(qPCRdata.raw, norm = "deltaCt", deltaCt.genes = rownames(qPCRdata.raw)[85:89], verbose = T)
plotCtDensity(dCT.norm)
plotCtDensity(rawData)
plotCtDensity(sr.norm)
plotCtDensity(q.norm)
plotCtBoxes(sr.norm, stratify = "type")
plotCtBoxes(dCT.norm, stratify = "type")
plotCtBoxes(q.norm, stratify = "type")

plotCtScatter(sr.norm, cards = c(1,2), col = "type", diag = T)
plotCtScatter(sr.norm, cards = c(3,4), col = "type", diag = T)
plotCtScatter(dCT.norm, cards = c(1,2), col = "type", diag = T)
plotCtScatter(dCT.norm, cards = c(3,4), col = "type", diag = T)

plotCtPairs(sr.norm, col = "type", diag = TRUE)
plotCtPairs(dCT.norm, col = "type", diag = TRUE)

plotCVBoxes(qPCRdata.raw, stratify = "type")

plotCtHeatmap(qPCRdata.raw, gene.names = "", dist = "euclidean")
plotCtHeatmap(dCT.norm, gene.names = "", dist = "euclidean")
clusterCt(dCT.norm, type = "samples")
clusterCt(dCT.norm, type = "genes")
cluster.list <- clusterCt(dCT.norm, type = "genes", n.cluster = 4, cex = 0.3)
cluster.list[[4]]

plotCtPCA(qPCRdata.raw, features = F)
groups <- as.factor(c("Control", "Control", "H2AZ_KD", "H2AZ_KD"))
qDE.ttest <- ttestCtData(dCT.norm, groups = groups)




