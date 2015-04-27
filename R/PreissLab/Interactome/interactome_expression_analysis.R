#  interactome_expression_analysis.R
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
#  

# load libraries
library("biomaRt")
library("gdata")
library("GO.db")
library("GEOquery")
library("ade4")

mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
filters <- listFilters(mouse)
attribs <- listAttributes(mouse)

interactome <- read.xls("~/Dropbox/REM project-Sebastian/HL-1 interactome superset.xlsx", sheet = "Sheet1")
colnames(interactome)[c(1,2)] <- c("ensembl_gene_id", "gene_symbol")

# GEO expression data (previously downloaded)
load("~/Data/Preiss/Interactome/gpl6246.rda")

# Cell paper on cardiac fibroblast re-programming
load("~/Data/Preiss/Interactome/gse22292.rda") 
load("~/Data/Preiss/Interactome/gse49192.rda")

interactome.add_id <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "affy_mogene_1_0_st_v1", "description", "external_gene_name"),
                            filters = "ensembl_gene_id",
                            values = interactome$ensembl_gene_id,
                            mart = mouse)
# remove NAs in Affy ID column
interactome.add_id <- interactome.add_id[!is.na(interactome.add_id$affy_mogene_1_0_st_v1), ]



# check which probes are available on the array
i1 <- as.character(intersect(rownames(gse22292$GSE22292_series_matrix.txt.gz@assayData$exprs), unique(interactome.add_id$affy_mogene_1_0_st_v1)))
# expression data from Cell paper (cardiac fibroblast transformation)
gse22292.exprs <- gse22292$GSE22292_series_matrix.txt.gz@assayData$exprs
gse49192.exprs <- gse49192[[2]]@assayData$exprs


# selecting following samples
# Neonatal fibroblasts, re-programmed fibroblast 2wk, re-programmed fibroblast 4wk, neonatal cardiomyocyte
s1 <- rownames(gse22292[[1]]@phenoData@data[c(16,17,18,4,5,6,7,8,9,13,14,15),])

neon.fib <- rownames(gse22292[[1]]@phenoData@data[c(16,17,18),])
rep.2wk.fib <- rownames(gse22292[[1]]@phenoData@data[c(4,5,6),])
rep.fail.2wk.fib <- rownames(gse22292[[1]]@phenoData@data[c(1,2,3),])
rep.4wk.fib <- rownames(gse22292[[1]]@phenoData@data[c(7,8,9),])
rep.fail.4wk.fib <- rownames(gse22292[[1]]@phenoData@data[c(10,11,12),])
neon.cardiomyo <- rownames(gse22292[[1]]@phenoData@data[c(13,14,15),])

interactome.exprs <- gse22292.exprs[i1, c(neon.fib, rep.fail.4wk.fib, rep.fail.2wk.fib, rep.2wk.fib, rep.4wk.fib, neon.cardiomyo)]
all.exprs <- gse22292.exprs[, c(neon.fib, rep.fail.4wk.fib, rep.fail.2wk.fib, rep.2wk.fib, rep.4wk.fib, neon.cardiomyo)]
interactome.exprs.sd <- apply(interactome.exprs, 1, sd)

# summarizing expression data on gene level,
l1 <- lapply(unique(interactome.add_id$ensembl_gene_id), function(x) unique(interactome.add_id[which(interactome.add_id$ensembl_gene_id == x), "affy_mogene_1_0_st_v1"])) 
names(l1) <- unique(interactome.add_id$ensembl_gene_id)
l2 <- lapply(l1, function(x) {
  r1 <- as.character(unlist(x))
  i1 <- as.character(intersect(rownames(interactome.exprs), r1))
  print(i1)
  if(length(i1) == 1) {
    avg <- interactome.exprs[i1, ]
  }
  else if(length(i1) > 1) {
    avg <- apply(interactome.exprs[i1, ], 2, mean)
    }
  }
  )
df1 <- t(data.frame(l2[lapply(l2, length)>0]))
# and on sample level?
interactome.exprs <- df1
interactome.exprs.sd <- apply(interactome.exprs, 1, sd)

#linear models for each of the interactome genes with expression data
fact <- gse22292[[1]]@phenoData@data[c(neon.fib, rep.fail.4wk.fib, rep.fail.2wk.fib, rep.2wk.fib, rep.4wk.fib, neon.cardiomyo), "characteristics_ch1.1"]
fact <- factor(fact, levels(fact)[c(5,2,4,1,3,6)])

fit1 <- apply(interactome.exprs, 1, function(x) lm(x ~ fact))
fit1.pvals <- lapply(fit1, function(x) anova(x)$'Pr(>F)'[1])
fit1.pvals.fdr <- p.adjust(unlist(fit1.pvals), method = "fdr")

# for all probes
fit2 <- apply(all.exprs, 1, function(x) lm(x ~ fact))
fit2.pvals <- lapply(fit2, function(x) anova(x)$'Pr(>F)'[1])
fit2.pvals.fdr <- p.adjust(unlist(fit2.pvals), method = "fdr")


#plotting of the data------------

pdf("Cardiac_development_boxplots.pdf", paper = "a4r")
par(mfrow = c(3,1))
lapply(names(which(interactome.exprs.sd > 1)), function(x){
  boxplot(interactome.exprs[x,] ~ fact, horizontal = F, pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5), axes = F, frame = F)
  axis(side = 1, lwd = 2.75, at = c(1,2,3,4,5,6), labels = c("CF", "2wk\nGFP-", "4wk\nGFP-", "2wk\niCMs", "4wk\niCMs", "CM"), padj = 1)
  axis(side = 2, lwd = 2.75, labels = T, cex.axis = 1.5)
  desc <- unique(interactome.add_id[which(interactome.add_id$ensembl_gene_id == x), "description"])
  desc <- unlist(lapply(strsplit(desc, "\\["), function(x) x[1]))
  title(paste(x, "\n", desc, sep = ""))
})
dev.off()

#expression analysis CF/2wkiCMs------------------
eset1 <- gse[[1]][, sampleNames(gse[[1]]) %in% c(rep.4wk.fib, neon.fib)]
pData(eset1)$new <- rep(NA, 6)
pData(eset1)$new <- c(rep("iCM_4wk", 3), rep("CF",3))
fact <- as.factor(pData(eset1)$new)
design <- model.matrix(~0+fact)
colnames(design) <- c("CF", "iCMs_4wk_vs_CF")
fit <- lmFit(eset1, design)
fit <- eBayes(fit)
topTable(fit, coef = ("iCMs_4wk_vs_CF"), adjust = "fdr")

sd1 <- apply(eset1@assayData$exprs, 1, sd)
