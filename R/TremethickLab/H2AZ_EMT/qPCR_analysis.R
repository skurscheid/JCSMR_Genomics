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
 "OCLN" = "ENSCAFG00000007805",
 "SIP1" = "ENSCAFG00000013859" ,
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
qPCRGenesTab$type <- "GOI"
qPCRGenes.control <- c("B2M", "GUSB", "HPRT1", "GAPDH", "LDHAL6B")
qPCRGenesTab[qPCRGenes.control,]$type <- "CONTROL"

#-----------Analysis of TGFb-treatment experiment--------------------
qPCRData.TGFb_exp1 <- read.csv("../../MDCK qPCR data/TGFb-treatment/20120620 PCR - TGFb-treatment Experiment 1.csv", header = T, row.names = 1, as.is = T)
qPCRData.TGFb_exp2 <- read.csv("../../MDCK qPCR data/TGFb-treatment/20120625 PCR - TGFb-treatment Experiment 2.csv", header = T, row.names = 1, as.is = T)
qPCRData.TGFb <- cbind(qPCRData.TGFb_exp1[, c(1:2)], qPCRData.TGFb_exp2[,c(1:2)], qPCRData.TGFb_exp1[, c(3:4)], qPCRData.TGFb_exp2[,c(3:4)])
colnames(qPCRData.TGFb) <- c(paste("Control_", c(1:4), sep= ""), paste("TGFb_treated_", c(1:4), sep= ""))

lapply(seq_along(1:ncol(qPCRData.TGFb)), function(x) {
  Ct <- qPCRData.TGFb[,x]
  genes <- rownames(qPCRData.TGFb)
  pos <- c(1:nrow(qPCRData.TGFb))
  flag <- rep("OK", nrow(qPCRData.TGFb))
  type <- c(rep("GOI", 84), rep("HKG", 5))
  fn <- paste("../../MDCK qPCR data/TGFb-treatment/", colnames(qPCRData.TGFb)[x], ".csv", sep = "")
  print(fn)
  tab <- data.frame(cbind(genes, Ct, pos, flag, type ))
  write.table(tab, file = fn, row.names = F, col.names = F, sep = "\t")
})

files <- data.frame(ID = 1:ncol(qPCRData.TGFb))
files$Files <- unlist(lapply(seq_along(1:ncol(qPCRData.TGFb)), function(x) { paste(colnames(qPCRData.TGFb)[x], ".csv", sep = "") }))

column.info <- list(flag = "flag", Ct = "Ct", feature = "genes", position = "pos", type = "type")

qPCRdata.raw <- readCtData(files$Files, 
                           format = "plain", 
                           path = "../../MDCK qPCR data/TGFb-treatment/", 
                           n.features = 89, 
                           column.info = list(flag = 4, Ct = 2, feature = 1, position = 3, type = 5), 
                           n.data = 1,
                           header = F)

pData(qPCRdata.raw) <- data.frame(Sample = c("Control", "Control", "Control", "Control", "TGFb-treated", "TGFb-treated", "TGFb-treated", "TGFb-treated"), Replicate = rep(1:4,2 ))

g <- featureNames(qPCRdata.raw)[1:10]
plotCtOverview(qPCRdata.raw, calibrator = "Control", genes = g, conf.int = T, ylim = c(0, 2), groups = c("Control", "Control", "Control", "Control", "TGFb-treated", "TGFb-treated", "TGFb-treated", "TGFb-treated"))

sr.norm <- normalizeCtData(qPCRdata.raw, norm = "scale.rank")
q.norm <- normalizeCtData(qPCRdata.raw, norm = "quantile")
dCT.norm <- normalizeCtData(qPCRdata.raw, norm = "deltaCt", deltaCt.genes = rownames(qPCRdata.raw)[85:89], verbose = T)

plotCtDensity(dCT.norm)
plotCtDensity(qPCRdata.raw)
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
cluster.list <- clusterCt(dCT.norm, type = "genes", n.cluster = 4, cex = 0.2)
cluster.list[[3]]

plotCtPCA(qPCRdata.raw, features = F)
groups <- as.factor(c("Control", "Control", "Control", "Control", "TGFb-treated", "TGFb-treated", "TGFb-treated", "TGFb-treated"))
qDE.ttest <- ttestCtData(dCT.norm, groups = groups)
qDE.ttest$genes <- as(qDE.ttest$genes, "character")
qDE.ttest$ensembl_gene_id <- qPCRGenesTab[qDE.ttest$genes, "ensembl_gene_id"]
qDE.ttest <- qDE.ttest[- which(is.na(qDE.ttest$ensembl_gene_id)), ]
rownames(qDE.ttest) <- qDE.ttest$ensembl_gene_id

#load("danpos2.anno.rda")
danpos2.anno[which(danpos2.anno$feature %in% qPCRGenesTab$ensembl_gene_id)]

#-----------Analysis of H2AZ RNAi experiment--------------------
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

files <- data.frame(ID = 1:ncol(qPCRData.H2AZ_KD))
files$Files <- unlist(lapply(seq_along(1:ncol(qPCRData.H2AZ_KD)), function(x) { paste(colnames(qPCRData.H2AZ_KD)[x], ".csv", sep = "") }))

column.info <- list(flag = "flag", Ct = "Ct", feature = "genes", position = "pos", type = "type")

qPCRdata.raw <- readCtData(files$Files, 
                           format = "plain", 
                           path = "../../MDCK qPCR data", 
                           n.features = 89, 
                           column.info = list(flag = 4, Ct = 2, feature = 1, position = 3, type = 5), 
                           n.data = 1,
                           header = F)

pData(qPCRdata.raw) <- data.frame(Sample = c("Control", "Control", "H2AZ_KD", "H2AZ_KD"), Replicate = rep(1:2,2 ))

g <- featureNames(qPCRdata.raw)[1:10]
plotCtOverview(qPCRdata.raw, calibrator = "Control", genes = g, conf.int = T, ylim = c(0, 2), groups = c("Control", "Control", "H2AZ_KD", "H2AZ_KD"))

sr.norm <- normalizeCtData(qPCRdata.raw, norm = "scale.rank")
q.norm <- normalizeCtData(qPCRdata.raw, norm = "quantile")
dCT.norm <- normalizeCtData(qPCRdata.raw, norm = "deltaCt", deltaCt.genes = rownames(qPCRdata.raw)[85:89], verbose = T)
plotCtDensity(dCT.norm)
plotCtDensity(qPCRdata.raw)
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
cluster.list <- clusterCt(dCT.norm, type = "genes", n.cluster = 4, cex = 0.2)
cluster.list[[3]]

plotCtPCA(qPCRdata.raw, features = F)
groups <- as.factor(c("Control", "Control", "H2AZ_KD", "H2AZ_KD"))
qDE.ttest <- ttestCtData(dCT.norm, groups = groups)




