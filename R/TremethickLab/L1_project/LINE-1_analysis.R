#  LINE-1_analysis.R
#  Copyright 2015 Sebastian Kurscheid <sebastian.kurscheid@anu.edu.au>
#
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
#------load libraries------
library("biomaRt")
library("GO.db")
library("GenomicFeatures")
library("GenomicAlignments")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("rtracklayer")
library("AnnotationDbi")
library("BSgenome.Hsapiens.UCSC.hg19")
library("metagene")
library("GenomicRanges")
library("RSQLite")
library("biovizBase")
library("Rsamtools")

#------environment------
host <- system("hostname", intern = T)
if (host %in% c("mmlab4.anu.edu.au")) {
  root <- "/Volumes/gduserv"
} else {
  root <- "/home/skurscheid"
}

setwd(paste(root, "/Data/Tremethick/LINE_1_project/", sep = ""))

#------some variables------
canonicalChr <- c(seq(1,19), "X", "Y", "M")

#------load annotation data------
repeatMasker <- read.table("~/Data/Annotations/GRCm38_mm10/repeatMaskerMM10.txt", header = T, as.is = T, sep = "\t")
repeatMasker$genoName <- gsub("chr", "", repeatMasker$genoName)
gr.repeatMasker <- GRanges(seqnames = repeatMasker$genoName, IRanges(start = repeatMasker$genoStart, end = repeatMasker$genoEnd), strand = repeatMasker$strand, mcols = repeatMasker[,c("repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft")])

#------subset LINE1 elements
gr.repeatMasker.LINE1 <- gr.repeatMasker[which(gr.repeatMasker$mcols.repFamily == "L1")]
seqlevels(gr.repeatMasker.LINE1, force = T) <- canonicalChr
gr.repeatMasker.LINE1 <- sort(gr.repeatMasker.LINE1)

#------load Ensembl TxDB------
mm10Ensembl <- loadDb(paste(root, "/Data/Annotations/GRCm38_mm10/mmusculus_gene_ensembl_GRCm38_TxDB.sqlite", sep = ""))
#------extract gene annotation
gr.genes <- genes(mm10Ensembl)
seqlevels(gr.genes, force = TRUE) <- canonicalChr
gr.genes <- sort(gr.genes)

#------import BAM files for visualization and quantitative analysis------------
# run picard MarkDuplicates prior to import to remove PCR duplicates
ga.dt_input <- import("DT_Input.sorted.rmdup.bam")
ga.dt_h2az <- import("DT_H2A.Z.sorted.rmdup.bam")
ga.g1_input <- import("G1_Input.sorted.rmdup.bam")
ga.g1_h2az <- import("G1_H2A.Z.sorted.rmdup.bam")
ga.s_input <- import("S_Input.sorted.rmdup.bam")
ga.s_h2az <- import("S_H2AZ.sorted.rmdup.bam")


seqlevels(ga.dt_input, force = TRUE) <- canonicalChr
seqlevels(ga.dt_h2az, force = TRUE) <- canonicalChr
seqlevels(ga.g1_input, force = TRUE) <- canonicalChr
seqlevels(ga.g1_h2az, force = TRUE) <- canonicalChr
seqlevels(ga.s_input, force = TRUE) <- canonicalChr
seqlevels(ga.s_h2az, force = TRUE) <- canonicalChr


gr.dt_input <- granges(ga.dt_input)
gr.dt_h2az <- granges(ga.dt_h2az)
gr.g1_input <- granges(ga.g1_input)
gr.g1_h2az <- granges(ga.g1_h2az)
gr.s_input <- granges(ga.s_input)
gr.s_h2az <- granges(ga.s_h2az)

width(gr.dt_input) <- 150
width(gr.dt_h2az) <- 150
width(gr.g1_input) <- 150
width(gr.g1_h2az) <- 150
width(gr.s_input) <- 150
width(gr.s_h2az) <- 150

cov.dt_input <- coverage(gr.dt_input)
cov.dt_h2az <- coverage(gr.dt_h2az)
cov.g1_input <- coverage(gr.g1_input)
cov.g1_h2az <- coverage(gr.g1_h2az)
cov.s_input <- coverage(gr.s_input)
cov.s_h2az <- coverage(gr.s_h2az)

#------results of using the peakranger "ranger" function (for narrow peaks)---------------------------------------------------
dat.peakranger_narrow <- data.frame(wd = rep("/home/skurscheid/Data/Tremethick/LINE_1_project/peakranger_analysis/", 6), 
                                    cell = c("DT", "DT", "G1", "G1", "S", "S"), 
                                    mark = rep(c("H2AZ", "H3K4me3"), 3), 
                                    peak = rep("narrow", 6))

grl.peakranger_narrow <- apply(dat.peakranger_narrow, 1, function(x) {
                               f <- paste(x["wd"], "/", 
                                           x["cell"], "_", 
                                           x["mark"], "_", 
                                           x["peak"], "/", 
                                           x["cell"], "_", 
                                           x["mark"], "_", 
                                           x["peak"], 
                                           "_region.bed", 
                                           sep = "")
                               if (file.exists(f) == TRUE){
                                gr <- import(f)
                                seqlevels(gr, force = TRUE) <- canonicalChr
                                }
                               return(gr)
                               })


grl.peakranger_narrow <- lapply(grl.peakranger_narrow, unstrand)
grl.peakranger_narrow <- lapply(grl.peakranger_narrow, sort)
grl.peakranger_narrow <- GRangesList(grl.peakranger_narrow)
names(grl.peakranger_narrow) <- apply(dat.peakranger_narrow, 1, function(x) paste(x["cell"], x["mark"], sep = "_"))

#------results of "arem" analysis (for narrow peaks)---------------------------------------------------
dat.arem_narrow <- data.frame(wd = rep("/home/skurscheid/Data/Tremethick/LINE_1_project/arem_analysis/", 3), 
                                    cell = c("DT", "G1", "S"), 
                                    mark = rep(c("H2AZ"), 3), 
                                    peak = rep("peaks", 3))

# use the XLS file as input, as statistics are only available here (FC, p-values, etc)
grl.arem_narrow <- apply(dat.arem_narrow, 1, function(x) {
  f <- paste(x["wd"], "/", 
             x["cell"], "_", 
             x["mark"], "/", 
             x["cell"], "_", 
             x["mark"], "_", 
             "arem", "_",
             x["peak"], 
             ".xls", 
             sep = "")
  if (file.exists(f) == TRUE){
    tab1 <- read.table(f, skip = 23, header = T, as.is = T, sep = "\t")
    gr <- GRanges(seqnames = tab1$chr, IRanges(start = tab1$start, end = tab1$end), strand = "*", mcols = tab1[, c("summit", "tags", "X.10.log10.pvalue.", "fold_enrichment", "FDR...")])
    seqlevels(gr, force = TRUE) <- canonicalChr
  }
  return(gr)
})

grl.arem_narrow <- lapply(grl.arem_narrow, unstrand)
grl.arem_narrow <- lapply(grl.arem_narrow, sort)
grl.arem_narrow <- GRangesList(grl.arem_narrow)
names(grl.arem_narrow) <- apply(dat.arem_narrow, 1, function(x) paste(x["cell"], x["mark"], sep = "_"))

# counting peaks intersecting with different L1 elements
dat.arem_narrow$L1_peaks <- rep(0,3)
dat.arem_narrow$L1Md_A_peaks <- rep(0,3)
dat.arem_narrow$L1_peaks <- sapply(names(grl.arem_narrow), function(x) length(subsetByOverlaps(grl.arem_narrow[[x]], gr.repeatMasker.LINE1)))
dat.arem_narrow$L1Md_A_peaks <- sapply(names(grl.arem_narrow), function(x) length(subsetByOverlaps(grl.arem_narrow[[x]], gr.repeatMasker.LINE1[grep("Md", gr.repeatMasker.LINE1$mcols.repName)])))

pdf("H2AZ_LINE1_L1Md_A_peak_counts.pdf", paper = "a4r")
ggplot(dat.arem_narrow, aes(cell, L1Md_A_peaks)) + geom_bar(stat = "identity") + ggtitle("H2AZ peaks at L1Md_A elements")
dev.off()

pdf("H2AZ_LINE1_L1_peak_counts.pdf", paper = "a4r")
ggplot(dat.arem_narrow, aes(cell, L1_peaks)) + geom_bar(stat = "identity") + ggtitle("H2AZ peaks at L1 elements")
dev.off()

# make histograms of fold-enrichment for all three groups - here all L1 elements
df1 <- data.frame(FE = subsetByOverlaps(grl.arem_narrow[[1]], gr.repeatMasker.LINE1)$mcols.fold_enrichment, source = rep("DT_H2AZ", length(subsetByOverlaps(grl.arem_narrow[[1]], gr.repeatMasker.LINE1))))
df1 <- rbind(df1, data.frame(FE = subsetByOverlaps(grl.arem_narrow[[2]], gr.repeatMasker.LINE1)$mcols.fold_enrichment, source = rep("G1_H2AZ", length(subsetByOverlaps(grl.arem_narrow[[2]], gr.repeatMasker.LINE1)))))
df1 <- rbind(df1, data.frame(FE = subsetByOverlaps(grl.arem_narrow[[3]], gr.repeatMasker.LINE1)$mcols.fold_enrichment, source = rep("S_H2AZ", length(subsetByOverlaps(grl.arem_narrow[[3]], gr.repeatMasker.LINE1)))))
df1$source <- factor(df1$source, levels(df1$source)[c(2,3,1)])

# histogram of fold enrichment at LINE1 element sites overlapping with H2AZ
pdf("H2AZ_L1_FE_histogram.pdf", paper = "a4r")
ggplot(df1, aes(x=FE, fill=source)) + geom_density(alpha=.3)+ ggtitle("Histogram of fold-enrichment at L1 associated H2AZ peaks")
dev.off()

# make histograms of fold-enrichment for all three groups - here only L1Md_A elements
gr1 <- gr.repeatMasker.LINE1[grep("Md", gr.repeatMasker.LINE1$mcols.repName)]
df2 <- data.frame(FE = subsetByOverlaps(grl.arem_narrow[[1]], gr1)$mcols.fold_enrichment, source = rep("DT_H2AZ", length(subsetByOverlaps(grl.arem_narrow[[1]], gr1))))
df2 <- rbind(df2, data.frame(FE = subsetByOverlaps(grl.arem_narrow[[2]], gr1)$mcols.fold_enrichment, source = rep("G1_H2AZ", length(subsetByOverlaps(grl.arem_narrow[[2]], gr1)))))
df2 <- rbind(df2, data.frame(FE = subsetByOverlaps(grl.arem_narrow[[3]], gr1)$mcols.fold_enrichment, source = rep("S_H2AZ", length(subsetByOverlaps(grl.arem_narrow[[3]], gr1)))))
df2$source <- factor(df2$source, levels(df2$source)[c(2,3,1)])

# histogram of fold enrichment at L1Md_A element sites overlapping with H2AZ
pdf("H2AZ_L1Md_A_FE_histogram.pdf", paper = "a4r")
ggplot(df2, aes(x=FE, fill=source)) + geom_density(alpha=.3) + ggtitle("Histogram of fold-enrichment at L1Md_A associated H2AZ peaks")
dev.off()


#------comparison between MACS2 and PeakRanger--------------------------------------------------------------------
# 
# G1 cells, H2A.Z, S cells
dat.macs2 <- data.frame(wd = as.character(rep("/home/skurscheid/Data/Tremethick/LINE_1_project/macs2_analysis/", 6)),
                        cell = c("DT", "DT", "G1", "G1", "S", "S"),
                        mark = rep(c("H2AZ", "H3K4me3"), 3),
                        peak = rep("narrow", 6),
                        suffix = rep("_peaks.narrowPeak", 6))

grl.macs2_narrow <- apply(dat.macs2, 1, function(x) {
                          f <- paste(x["wd"], "/", x["cell"], "_", x["mark"], "_", x["peak"], "/", x["cell"], "_", x["mark"], "_", x["peak"], x["suffix"], sep = "")
                          print(f)
                          if (file.exists(f)){
                            t1 <- read.table(f, as.is = T, header = F, sep = "\t")
                            gr <- GRanges(seqnames = t1$V1, IRanges(start = t1$V2, end = t1$V3), strand = "*", mcols = t1[,c("V4", "V5", "V7", "V8", "V9", "V10")])
                            seqlevels(gr, force = TRUE) <- canonicalChr
                          }
                          return(gr)
                        })

grl.macs2_narrow <- lapply(grl.macs2_narrow, unstrand)
grl.macs2_narrow <- lapply(grl.macs2_narrow, sort)
grl.macs2_narrow <- GRangesList(grl.macs2_narrow)
names(grl.macs2_narrow) <- apply(dat.macs2, 1, function(x) paste(x["cell"], x["mark"], sep = "_"))



# overlap between MACS2 and PeakRanger results
grl.consensus <- lapply(seq_along(1:nrow(dat.macs2)), function(x) {
                        gr1 <- unlist(grl.peakranger_narrow[x])
                        gr2 <- unlist(grl.macs2_narrow[x])
                        subsetByOverlaps(gr1, gr2)})

names(grl.consensus) <- apply(dat.macs2, 1, function(x) paste(x["cell"], x["mark"], sep = "_"))
grl.consensus <- GRangesList(grl.consensus)


# check which peaks overlap with LINE1 elements
grl.consensus.LINE1 <- lapply(grl.consensus, function(x) {
                             gr <- subsetByOverlaps(gr.repeatMasker.LINE1, unlist(x))
                             return(gr)
                             })

grl.consensus.LINE1 <- GRangesList(grl.consensus.LINE1)

grl.consensus.h2az.LINE1 <- lapply(grl.consensus[c(1,3,5)], function(x) {
                                   gr <- subsetByOverlaps(unlist(x), gr.repeatMasker.LINE1)
                                   return(gr)
                                  })
grl.consensus.h2az.LINE1 <- GRangesList(grl.consensus.h2az.LINE1)
grl.consensus.h2az.LINE1 <- resize(grl.consensus.h2az.LINE1, 1000, fix = "center")
grl.consensus.h2az.LINE1 <- GRangesList(lapply(grl.consensus.h2az.LINE1, trim))
grl.consensus.h2az.LINE1 <- GRangesList(lapply(grl.consensus.h2az.LINE1, sort))

#------DT data----------------------------------------------

grl.dt.h2az.LINE1 <- GRangesList(lapply(canonicalChr, function(x) grl.consensus.h2az.LINE1$"DT_H2AZ"[which(seqnames(grl.consensus.h2az.LINE1$"DT_H2AZ") == x)]))
irl.dt.h2az.LINE1 <- as(grl.dt.h2az.LINE1, "IRangesList")
vw.dt.h2az.LINE1 <- Views(cov.dt_h2az, irl.dt.h2az.LINE1)
vw.dt.input.LINE1 <- Views(cov.dt_input, irl.dt.h2az.LINE1)

export(flatGrl(grl.dt.h2az.LINE1), format = "BED", con = "dt.h2az.line1.bed")

l1 <- lapply(names(vw.dt.h2az.LINE1), function(x) {
  e1 <- unlist(vw.dt.h2az.LINE1[x])
  if (length(e1) > 0){
    df1 <- data.frame(lapply(e1, function(y) as(unlist(y), "integer")))
    if (nrow(df1) == 1000){
      return(df1)  
    }
  }
  }
)
l1 <- l1[-(which(sapply(l1,is.null),arr.ind=TRUE))]
df.dt.h2az.LINE1 <- data.frame(l1)

l1 <- lapply(names(vw.dt.input.LINE1), function(x) {
  e1 <- unlist(vw.dt.input.LINE1[x])
  if (length(e1) > 0){
    df1 <- data.frame(lapply(e1, function(y) as(unlist(y), "integer")))
    if (nrow(df1) == 1000){
      return(df1)  
    }
  }
}
)
l1 <- l1[-(which(sapply(l1,is.null),arr.ind=TRUE))]
df.dt.input.LINE1 <- data.frame(l1)

pdf("DT_H2AZ.pdf")
plot(-499:500, apply(df.dt.h2az.LINE1[1:225], 1, mean), type = "l")
abline(v = lowess(apply(df.dt.h2az.LINE1, 1, mean) ~ -499:500)$x[which.max(lowess(apply(df.dt.h2az.LINE1, 1, mean) ~ -499:500)$y)], col = "grey", lty = 2)
lines(lowess(apply(df.dt.h2az.LINE1[1:225], 1, mean) ~ -499:500),col="green3")
dev.off()

pdf("DT_Input.pdf")
plot(-499:500, apply(df.dt.input.LINE1, 1, mean), type = "l")
lines(lowess(apply(df.dt.input.LINE1[1:225], 1, mean) ~ -499:500),col="green3")
dev.off()

#------G1 data----------------------------------------------
grl.g1.h2az.LINE1 <- GRangesList(lapply(canonicalChr, function(x) grl.consensus.h2az.LINE1$"G1_H2AZ"[which(seqnames(grl.consensus.h2az.LINE1$"G1_H2AZ") == x)]))
irl.g1.h2az.LINE1 <- as(grl.g1.h2az.LINE1, "IRangesList")
vw.g1.h2az.LINE1 <- Views(cov.g1_h2az, irl.g1.h2az.LINE1)
vw.g1.input.LINE1 <- Views(cov.g1_input, irl.g1.h2az.LINE1)

export(flatGrl(grl.g1.h2az.LINE1), format = "BED", con = "g1.h2az.line1.bed")

l1 <- lapply(names(vw.g1.h2az.LINE1), function(x) {
  e1 <- unlist(vw.g1.h2az.LINE1[x])
  if (length(e1) > 0){
    df1 <- data.frame(lapply(e1, function(y) as(unlist(y), "integer")))
    if (nrow(df1) == 1000){
      return(df1)  
    }
  }
}
)
l1 <- l1[-(which(sapply(l1,is.null),arr.ind=TRUE))]
df.g1.h2az.LINE1 <- data.frame(l1)

l1 <- lapply(names(vw.g1.input.LINE1), function(x) {
  e1 <- unlist(vw.g1.input.LINE1[x])
  if (length(e1) > 0){
    df1 <- data.frame(lapply(e1, function(y) as(unlist(y), "integer")))
    if (nrow(df1) == 1000){
      return(df1)  
    }
  }
}
)
l1 <- l1[-(which(sapply(l1,is.null),arr.ind=TRUE))]
df.g1.input.LINE1 <- data.frame(l1)

pdf("G1_H2AZ.pdf")
plot(-499:500, apply(df.g1.h2az.LINE1, 1, mean), type = "l")
abline(v = lowess(apply(df.g1.h2az.LINE1, 1, mean) ~ -499:500)$x[which.max(lowess(apply(df.g1.h2az.LINE1, 1, mean) ~ -499:500)$y)], col = "grey", lty = 2)
lines(lowess(apply(df.g1.h2az.LINE1, 1, mean) ~ -499:500),col="green3")
dev.off()

pdf("G1_Input.pdf")
plot(-499:500, apply(df.g1.input.LINE1, 1, mean), type = "l")
lines(lowess(apply(df.g1.input.LINE1, 1, mean) ~ -499:500),col="green3")
dev.off()

#------for S phase data----------------------------------------------
grl.s.h2az.LINE1 <- GRangesList(lapply(canonicalChr, function(x) grl.consensus.h2az.LINE1$"S_H2AZ"[which(seqnames(grl.consensus.h2az.LINE1$"S_H2AZ") == x)]))
irl.s.h2az.LINE1 <- as(grl.s.h2az.LINE1, "IRangesList")
vw.s.h2az.LINE1 <- Views(cov.s_h2az, irl.s.h2az.LINE1)
vw.s.input.LINE1 <- Views(cov.s_input, irl.s.h2az.LINE1)

# Annotation of LINE-1 elements - crude methods
gr1 <- unlist(grl.s.h2az.LINE1)
l1 <- lapply(gr1, function(x) subsetByOverlaps(gr.genes, x))
gr1$gene_id <- as.character(lapply(l1, function(x) x$gene_id))
l1 <- lapply(gr1, function(x) select(mm10, keys=x$gene_id, keytype="ensembl_gene_id", columns=c("description")))
gr1$description <- as.character(l1)


export(flatGrl(grl.s.h2az.LINE1), format = "BED", con = "s.h2az.line1.bed")

l1 <- lapply(names(vw.s.h2az.LINE1), function(x) {
  e1 <- unlist(vw.s.h2az.LINE1[x])
  if (length(e1) > 0){
    df1 <- data.frame(lapply(e1, function(y) as(unlist(y), "integer")))
    if (nrow(df1) == 1000){
      return(df1)  
    }
  }
}
)

l1 <- l1[-(which(sapply(l1,is.null),arr.ind=TRUE))]
df.s.h2az.LINE1 <- data.frame(l1)

l1 <- lapply(names(vw.s.input.LINE1), function(x) {
  e1 <- unlist(vw.s.input.LINE1[x])
  if (length(e1) > 0){
    df1 <- data.frame(lapply(e1, function(y) as(unlist(y), "integer")))
    if (nrow(df1) == 1000){
      return(df1)  
    }
  }
}
)

l1 <- l1[-(which(sapply(l1,is.null),arr.ind=TRUE))]
df.s.input.LINE1 <- data.frame(l1)

pdf("S_H2AZ.pdf")
plot(-499:500, apply(df.s.h2az.LINE1, 1, mean), type = "l")
abline(v = lowess(apply(df.s.h2az.LINE1, 1, mean) ~ -499:500)$x[which.max(lowess(apply(df.s.h2az.LINE1, 1, mean) ~ -499:500)$y)], col = "grey", lty = 2)
lines(lowess(apply(df.s.h2az.LINE1, 1, mean) ~ -499:500),col="green3")
dev.off()

pdf("S_Input.pdf")
plot(-499:500, apply(df.s.input.LINE1, 1, mean), type = "l")
lines(lowess(apply(df.s.input.LINE1, 1, mean) ~ -499:500),col="green3")
dev.off()


#---------------------------------------------------------------------------
grl.consensus.h3k4me3.LINE1 <- lapply(grl.consensus[c(2,4)], function(x) {
                                      gr <- subsetByOverlaps(unlist(x), gr.repeatMasker.LINE1)
                                      return(gr)
                                      })
grl.consensus.h3k4me3.LINE1 <- GRangesList(grl.consensus.h3k4me3.LINE1)
grl.consensus.h3k4me3.LINE1 <- resize(grl.consensus.h3k4me3.LINE1, 1000, fix = "center")
irl.consensus.h3k4me3.LINE1 <- as(grl.consensus.h3k4me3.LINE1, "IRangesList")



gr.common.h2az.LINE1 <- subsetByOverlaps(unlist(grl.consensus.LINE1[1]), unlist(grl.consensus.LINE1[3]))

grl.consensus.LINE1.flank1kb <- lapply(grl.consensus, function(x) {
  gr <- subsetByOverlaps(flank(gr.repeatMasker.LINE1, 1000, both = TRUE), unlist(x))
  return(gr)
})

grl.consensus.LINE1.flank1kb <- GRangesList(grl.consensus.LINE1.flank1kb)

gr.common.h2az.LINE1.flank1kb <- subsetByOverlaps(unlist(grl.consensus.LINE1.flank1kb[1]), unlist(grl.consensus.LINE1.flank1kb[3]))
gr.common.h3k4me3.LINE1.flank1kb <- subsetByOverlaps(unlist(grl.consensus.LINE1.flank1kb[2]), unlist(grl.consensus.LINE1.flank1kb[4]))

# extract LINE1 associated H2A.Z peaks which are unique to each cell line
# perhaps a bit crude, but I can't think of a better option right now
# essentially I am extracting the gaps from the regions common to both cell lines and 
# overlap them with the LINE1 peaks in one cell line
# this misses those chromosomes where there are no common peaks, 
# so a bit of post-processing is required

# DT
gr.consensus.LINE1.dt_h2az <- subsetByOverlaps(grl.consensus.LINE1.flank1kb[[1]], gaps(gr.common.h2az.LINE1.flank1kb))
table(seqnames(grl.consensus.LINE1.flank1kb[[1]]))
#  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19   X   Y   M 
#135 174  63  71  94  75 108  53  50  38 162  34  75  34  35  40  72  33  38   7   0   0 
table(seqnames(gr.common.h2az.LINE1.flank1kb))
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  X  Y  M 
# 41 45  4  8 22 14 19 14  8  6 30  5 20  1  8  8 18  8  1  0  0  0 
# so peaks on X are missing
gr.consensus.LINE1.dt_h2az <- c(gr.consensus.LINE1.dt_h2az, grl.consensus.LINE1.flank1kb[[1]][which(seqnames(grl.consensus.LINE1.flank1kb[[1]]) %in% c("X"))])
table(seqnames(gr.consensus.LINE1.dt_h2az))
gr.consensus.LINE1.dt_h2az <- sort(gr.consensus.LINE1.dt_h2az)
grl.consensus.LINE1.dt_h2az <- GRangesList(lapply(canonicalChr, function(x) gr.consensus.LINE1.dt_h2az[which(seqnames(gr.consensus.LINE1.dt_h2az) == x)]))
irl.consensus.LINE1.dt_h2az <- as(grl.consensus.LINE1.dt_h2az, "IRangesList")
# G1
gr.consensus.LINE1.g1_h2az <- subsetByOverlaps(grl.consensus.LINE1.flank1kb[[3]], gaps(gr.common.h2az.LINE1.flank1kb))
table(seqnames(grl.consensus.LINE1.flank1kb[[3]]))
#  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19   X   Y   M 
#106  92  28  39  60  52  52  53  49  39  71  26  68  25  23  25  31  23  20   4   0   0 

table(seqnames(gr.common.h2az.LINE1.flank1kb))
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  X  Y  M 
# 44 45  7  9 26 18 17 15 11 10 33  5 23  0  7  9 16  8  2  0  0  0 
# so peaks on X are missing
gr.consensus.LINE1.g1_h2az <- c(gr.consensus.LINE1.g1_h2az, grl.consensus.LINE1.flank1kb[[3]][which(seqnames(grl.consensus.LINE1.flank1kb[[3]]) %in% c("X"))])
table(seqnames(gr.consensus.LINE1.g1_h2az))
gr.consensus.LINE1.g1_h2az <- sort(gr.consensus.LINE1.g1_h2az)
grl.consensus.LINE1.g1_h2az <- GRangesList(lapply(canonicalChr, function(x) gr.consensus.LINE1.g1_h2az[which(seqnames(gr.consensus.LINE1.g1_h2az) == x)]))
irl.consensus.LINE1.g1_h2az <- as(grl.consensus.LINE1.g1_h2az, "IRangesList")

# important to keep in mind that "consensus" here refers to LINE1 elements which had associated peaks
# identified by both MACS2 and peakranger in the respective sample
vw.dt_h2az <- Views(cov.dt_h2az, irl.consensus.LINE1.dt_h2az)
vw.dt_input <- Views(cov.dt_input, irl.consensus.LINE1.dt_h2az)
vw.g1_h2az <- Views(cov.g1_h2az, irl.consensus.LINE1.g1_h2az)
vw.g1_input <- Views(cov.g1_input, irl.consensus.LINE1.g1_h2az)




