#  tcga_esca_cimp_h2abbd_analysis.R
#
#  Copyright 2015 Sebastian Kurscheid <sebastian.kurscheid@anu.edu.au>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

source("/Users/u1001407/Dropbox/Development/GeneralPurpose/R/heatmap.3.R")

library(ade4)
library(cluster)
library(matlab)
library(GenomicRanges)

load("~/Data/Annotations/Platforms/gr.annot450k.rda")
gr1 <- gr.annot450k[which(gr.annot450k$Probe_SNPs == "" & gr.annot450k$Probe_SNPs_10 == "")]
seqlevels(gr1, force = T) <- as.character(seq(1,22,1))
gr.annot450k.auto.noSNPs <- gr1

#--------------------------------------------------------------------------------
#
# HM450 data
#
#--------------------------------------------------------------------------------

# load data 
esca.hm450.metadata <- read.table("METADATA/", header = T, as.is = T, sep = "\t")
esca.hm450.path <- ("DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/")
esca.hm450.files <- list.files("DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/")
esca.hm450.samples <- unlist(lapply(strsplit(unlist(lapply(strsplit(esca.hm450.files, "\\."), function(x) x[6])), "-"), function(y) paste(y[1], y[2], y[3], sep = "-")))

#-----------------------------------------------------------------
for (x in esca.hm450.files) {
  print(x)
  t1 <- read.table(paste(esca.hm450.path, x, sep = ""), header = T, as.is = T, sep = "\t")
  if (which(esca.hm450.files == x) == 1) {
    c1 <- t1[,1]
    beta1 <- t1[,2]
    t2 <- cbind(c1,beta1)
  } else {
    t2 <- cbind(t2, t1[,2])
  }
}

t2 <- t2[-1,]
rownames(t2) <- t2[,1]
t2 <- [,-1]
colnames(t2) <- esca.hm450.samples
esca.hm450 <- t2
rm(t2)

esca.hm450 <- esca.hm450[names(gr.annot450k.auto.noSNPs),]
esca.hm450 <- esca.hm450[complete.cases(esca.hm450),]

#--------------------------------------------------------------------------------
#
# RNA-Seq data
#
#--------------------------------------------------------------------------------

# load data
esca.rnaseq.metadata <- read.table("METADATA/BCGSC__IlluminaHiSeq_RNASeq/bcgsc.ca_ESCA.IlluminaHiSeq_RNASeq.sdrf.txt", header = T, as.is = T, sep = "\t")
esca.rnaseq.samples <- unlist(lapply(strsplit(esca.rnaseq.metadata[which(esca.rnaseq.metadata$Comment..TCGA.Data.Type..1 == "Expression-Gene"), "Comment..TCGA.Barcode."], "-"), function(x) paste(x[1], x[2], x[3], sep = "-")))
rnaseq.files <- esca.rnaseq.metadata[which(esca.rnaseq.metadata$Comment..TCGA.Data.Type..1 == "Expression-Gene"),]$Derived.Data.File
rnaseq.path <- ("RNASeq/BCGSC__IlluminaHiSeq_RNASeq/Level_3/")

# here we use the RSEM estimates
for (x in rnaseq.files) {
  print(x)
  f1 <- paste(rnaseq.path, x, sep = "")
  if (file.exists(f1)) {
    t1 <- read.table(f1, header = T, as.is = T, sep = "\t")
    if (which(rnaseq.files == x) == 1) {
      c1 <- t1[,1]
      scaled_estimates <- t1[,"RPKM"]
      t2 <- cbind(c1, scaled_estimates)
    } else {
      t2 <- cbind(t2, t1[, "RPKM"])
    }
  }
}

rownames(t2) <- t2[,1]
t2 <- t2[,-1]
class(t2) <- "numeric"

colnames(t2) <- esca.rnaseq.samples
esca.rnaseq <- t2
rm(t2)

esca.rnaseq <- esca.rnaseq[,-which(duplicated(colnames(esca.rnaseq)))]
esca.rnaseq.samples <- esca.rnaseq.samples[-which(duplicated(esca.rnaseq.samples))]

# here we use the RSEM normalized counts
for (x in rnaseq.files.rsem_norm) {
  print(x)
  f1 <- paste(rnaseq.path, x, sep = "")
  if (file.exists(f1)) {
    t1 <- read.table(f1, header = T, as.is = T, sep = "\t")
    if (which(rnaseq.files.rsem_norm == x) == 1) {
      c1 <- t1[,1]
      norm_counts <- t1[,"normalized_count"]
      t2 <- cbind(c1, norm_counts)
    } else {
      t2 <- cbind(t2, t1[, "normalized_count"])
    }
  }
}

rownames(t2) <- t2[,1]
t2 <- t2[,-1]
class(t2) <- "numeric"

esca.rnaseq.samples <- unlist(lapply(strsplit(esca.rnaseq.metadata[which(esca.rnaseq.metadata$Comment..TCGA.Data.Type..1 == "RSEM_genes"), "Comment..TCGA.Barcode."], "-"), function(x) paste(x[1], x[2], x[3], sep = "-")))

colnames(t2) <- esca.rnaseq.samples
esca.rnaseq.rsem_norm <- t2
rm(t2)
esca.rnaseq.rsem_norm <- esca.rnaseq.rsem_norm[,-which(duplicated(esca.rnaseq.samples))]
esca.rnaseq.samples <- esca.rnaseq.samples[-which(duplicated(esca.rnaseq.samples))]

#-----------------------------------------------------------------
# create tables with samples shared between both data sets
i1 <- intersect(esca.hm450.samples, esca.rnaseq.samples)

esca.rnaseq <- esca.rnaseq[,i1]
esca.hm450 <- esca.hm450[,i1]

# select top 1000 probes (based on SD)
sd1 <- apply(esca.hm450, 1, sd)
sd1.top1k <- names(sd1[order(sd1, decreasing = T)][1:1000])

# perform PCA for normalizing data
pca1 <- dudi.pca(t(esca.hm450[sd1.top1k,]), scale = TRUE, center = TRUE, scannf = F, nf = 5)
pdf("heatmap.450k.pdf")
hm1 <- heatmap.3(t(as.matrix(pca1$tab)), trace = "none", hclustfun=function(x) hclust(x,method="ward"), col=jet.colors(10))
dev.off()

# classify CIMP status
tab <- pca1$tab
hccol <- as.hclust(hm1$colDendrogram)
silhouette1 <- list()

for(i in 1:9)
  silhouette1[[i]] <- silhouette(cutree(hccol, k = i+1), daisy(tab))

pdf("shilhouetteHeatAutosome.pdf",pointsize=10, height=8,width=8)
par(mfrow=c(3,3))
for(i in 1:9)
  plot(silhouette1[[i]])
dev.off()

AutosomeCluster <- as.data.frame(cbind(autok2=cutree(hccol, k = 2),autok3=cutree(hccol, k = 3),autok4=cutree(hccol, k = 4)))
esca.cimp1 <- rownames(AutosomeCluster[which(AutosomeCluster[,1] == 1),])
esca.cimp2 <- rownames(AutosomeCluster[which(AutosomeCluster[,1] == 2),])

rc.cimp <- c("grey", "black")[match(AutosomeCluster[,1], c(1,2))]
rcw <- cbind(rc.cimp, rc.cimp)

rc.cimp2 <- jet.colors(3)[match(AutosomeCluster[,2], c(1,2,3))]
rcw <- cbind(rc.cimp, rc.cimp2)
esca.cimp1 <- rownames(AutosomeCluster[which(AutosomeCluster[,2] == 1),])
esca.cimp2 <- rownames(AutosomeCluster[which(AutosomeCluster[,2] == 2),])
esca.cimp3 <- rownames(AutosomeCluster[which(AutosomeCluster[,2] == 3),])

pdf("heatmap.esca.450k.cimp.pdf", height = 8, width = 8)
HeatAutosomeProbes <- heatmap.3(t(tab),col=jet.colors(10),
                                trace="none",density="density",denscol="black",
                                hclustfun=function(x) hclust(x,method="ward"),
                                ColSideColors=rcw,cex.var=0.2,
                                scale="none",xlab="sample",ylab="probe",
                                labRow = rep("",ncol(tab)),labCol = rep("",nrow(tab)),
                                main="Normalized DNA methylation (autosomal probes)")
legend("topright", legend = c("CIMP1", "CIMP2"), fill = c("grey", "black"), ncol = 4, cex = 0.5)
dev.off()


