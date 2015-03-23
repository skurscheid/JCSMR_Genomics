#  tcga_coad_cimp_h2abbd_analysis.R
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

library(ade4)
library(cluster)
library(matlab)
library(GenomicRanges)

load("~/Data/Annotations/Platforms/gr.annot450k.rda")
gr1 <- gr.annot450k[which(gr.annot450k$Probe_SNPs == "" & gr.annot450k$Probe_SNPs_10 == "")]
seqlevels(gr1, force = T) <- as.character(seq(1,22,1))
gr.annot450k.auto.noSNPs <- gr1

coad.rnaseq.metadata <- read.table("METADATA/UNC__IlluminaHiSeq_RNASeqV2/unc.edu_COAD.IlluminaHiSeq_RNASeqV2.1.12.0.sdrf.txt", header = T, as.is = T, sep = "\t")
coad.rnaseq.samples <- unlist(lapply(strsplit(coad.rnaseq.samples[which(coad.rnaseq.metadata$Comment..TCGA.Data.Type..1 == "RSEM_genes"), "Comment..TCGA.Barcode."], "-"), function(x) paste(x[1], x[2], x[3], sep = "-")))
rnaseq.files <- coad.rnaseq.metadata[which(coad.rnaseq.metadata$Comment..TCGA.Data.Type..1 == "RSEM_genes"),]$Derived.Data.File
rnaseq.files.norm <- coad.rnaseq.metadata[which(coad.rnaseq.metadata$Comment..TCGA.Data.Type..1 == "RSEM_genes_normalized"),]$Derived.Data.File
rnaseq.path <- ("RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/")

# here we use the RSEM estimates
for (x in rnaseq.files) {
  print(x)
  f1 <- paste(rnaseq.path, x, sep = "")
  if (file.exists(f1)) {
    t1 <- read.table(f1, header = T, as.is = T, sep = "\t")
    if (which(rnaseq.files == x) == 1) {
      c1 <- t1[,1]
      scaled_estimates <- t1[,3]
      t2 <- cbind(c1, scaled_estimates)
    } else {
      t2 <- cbind(t2, t1[, 3])
    }
  }
}

rownames(t2) <- t2[,1]
t2 <- [,-1]
class(t2) <- "numeric"

coad.rnaseq.samples <- coad.rnaseq.samples[-which(sapply(rnaseq.files.norm, function(x) file.exists(paste(rnaseq.path, x, sep = ""))) == FALSE)]
colnames(t2) <- coad.rnaseq.samples
coad.rnaseq <- t2

# here we use the normalized RSEM estimates
for (x in rnaseq.files.norm) {
  print(x)
  f1 <- paste(rnaseq.path, x, sep = "")
  if (file.exists(f1)) {
    t1 <- read.table(f1, header = T, as.is = T, sep = "\t")
    if (which(rnaseq.files.norm == x) == 1) {
      c1 <- t1[,1]
      norm_counts <- t1[,2]
      t2 <- cbind(c1, norm_counts)
    } else {
      t2 <- cbind(t2, t1[, 2])
    }
  }
}

rownames(t2) <- t2[,1]
t2 <- t2[,-1]
class(t2) <- "numeric"

coad.rnaseq.norm_count.samples <- unlist(lapply(strsplit(coad.rnaseq.metadata[which(coad.rnaseq.metadata$Comment..TCGA.Data.Type..1 == "RSEM_genes_normalized"), "Comment..TCGA.Barcode."], "-"), function(x) paste(x[1], x[2], x[3], sep = "-")))
coad.rnaseq.norm_count.samples <- coad.rnaseq.norm_count.samples[-which(sapply(rnaseq.files.norm, function(x) file.exists(paste(rnaseq.path, x, sep = ""))) == FALSE)]
colnames(t2) <- coad.rnaseq.samples
coad.rnaseq.norm_count <- t2

#-----------------------------------------------------------------
for (x in hm450.files) {
  print(x)
  t1 <- read.table(paste(hm450.path, x, sep = ""), header = T, as.is = T, sep = "\t")
  if (which(hm450.files == x) == 1) {
    c1 <- t1[,1]
    beta1 <- t1[,2]
    t2 <- cbind(c1,beta1)
  } else {
    t2 <- cbind(t2, t1[,2])
  }
}


coad.hm450.samples <- (unlist(lapply(strsplit(hm450.samples, "-"), function(x) paste(x[1], x[2], x[3], sep = "-"))))
colnames(t2) <- coad.hm450.samples
class(t2) <- "numeric"

i1 <- intersect(coad.hm450.samples, coad.rnaseq.samples)

coad.rnaseq <- coad.rnaseq[,i1]

coad.hm450 <- t2[names(gr.annot450k.auto.noSNPs),]
coad.hm450 <- coad.hm450[complete.cases(coad.hm450),]
coad.hm450 <- coad.hm450[,-which(duplicated(colnames(coad.hm450)))]
coad.hm450 <- coad.hm450[,i1]

# select top 1000 probes (based on SD)
sd1 <- apply(coad.hm450, 1, sd)
sd1.top1k <- names(sd1[order(sd1, decreasing = T)][1:1000])

# perform PCA for normalizing data
pca1 <- dudi.pca(t(coad.hm450[sd1.top1k,]), scale = TRUE, center = TRUE, scannf = F, nf = 5)
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
coad.cimp1 <- rownames(AutosomeCluster[which(AutosomeCluster[,1] == 1),])
coad.cimp2 <- rownames(AutosomeCluster[which(AutosomeCluster[,1] == 2),])

rc.cimp <- c("grey", "black")[match(AutosomeCluster[,1], c(1,2))]
rcw <- cbind(rc.cimp, rc.cimp)

pdf("heatmap.coad.450k.cimp.pdf", height = 8, width = 8)
HeatAutosomeProbes <- heatmap.3(t(tab),col=jet.colors(10),
                                trace="none",density="density",denscol="black",
                                hclustfun=function(x) hclust(x,method="ward"),
                                ColSideColors=rcw,cex.var=0.2,
                                scale="none",xlab="sample",ylab="probe",
                                labRow = rep("",ncol(tab)),labCol = rep("",nrow(tab)),
                                main="Normalized DNA methylation (autosomal probes)")
legend("topright", legend = c("CIMP1", "CIMP2"), fill = c("grey", "black"), ncol = 4, cex = 0.5)
dev.off()


