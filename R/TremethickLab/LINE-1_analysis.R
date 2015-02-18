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
# load libraries
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

# environment
setwd("/home/skurscheid/Data/Tremethick/LINE_1_project/")

# some variables
canonicalChr <- c(seq(1,19), "X", "Y", "M")

# load annotation data
repeatMasker <- read.table("~/Data/Annotations/GRCm38_mm10/repeatMaskerMM10.txt", header = T, as.is = T, sep = "\t")
gr.repeatMasker <- GRanges(seqnames = repeatMasker$genoName, IRanges(start = repeatMasker$genoStart, end = repeatMasker$genoEnd), strand = repeatMasker$strand, mcols = repeatMasker[,c("repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft")])
# subset LINE1 elements
gr.repeatMasker.LINE1 <- gr.repeatMasker[which(gr.repeatMasker$mcols.repClass == "LINE")]
seqlevels(gr.repeatMasker.LINE1, force = T) <- canonicalChr
gr.repeatMasker.LINE1 <- sort(gr.repeatMasker.LINE1)

# load peakranger ccat results - detecting broad regions
gr.peakranger_dt_h2az <- import("peakranger_analysis/DT_H2A.Z/DT_H2A.Z_region.bed")
gr.peakranger_dt_h3k4me3 <- import("peakranger_analysis/DT_H3K4me3/DT_H3K4me3_region.bed")
gr.peakranger_g1_h2az <- import("peakranger_analysis/G1_H2A.Z/G1_H2A.Z_region.bed")
gr.peakranger_g1_h3k4me3 <- import("peakranger_analysis/G1_H3K4me3/G1_H3K4me3_region.bed")

# only work with canonical chromosomes
seqlevels(gr.repeatMasker, force = T) <- canonicalChr
seqlevels(gr.peakranger_dt_h2az, force = T) <- canonicalChr
seqlevels(gr.peakranger_g1_h2az, force = T) <- canonicalChr
seqlevels(gr.peakranger_g1_h3k4me3, force = T) <- canonicalChr
seqlevels(gr.peakranger_dt_h3k4me3, force = T) <- canonicalChr

# only retain regions with FDR <= 10%
gr.peakranger_dt_h2az <- gr.peakranger_dt_h2az[which(gr.peakranger_dt_h2az$score <= 0.10)]
gr.peakranger_g1_h2az <- gr.peakranger_g1_h2az[which(gr.peakranger_g1_h2az$score <= 0.10)]
gr.peakranger_g1_h3k4me3 <- gr.peakranger_g1_h3k4me3[which(gr.peakranger_g1_h3k4me3$score <= 0.10)]
gr.peakranger_dt_h3k4me3 <- gr.peakranger_dt_h3k4me3[which(gr.peakranger_dt_h3k4me3$score <= 0.10)]

# inspection of data reveals that only "+" strand regions present
# this seems odd, for now I will replace "+" with "*"
strand(gr.peakranger_dt_h2az) <- "*"
strand(gr.peakranger_g1_h2az) <- "*"
strand(gr.peakranger_g1_h3k4me3) <- "*"
strand(gr.peakranger_dt_h3k4me3) <- "*"

# sort the data
gr.peakranger_dt_h2az <- sort(gr.peakranger_dt_h2az)
gr.peakranger_g1_h2az <- sort(gr.peakranger_g1_h2az)
gr.peakranger_g1_h3k4me3 <- sort(gr.peakranger_g1_h3k4me3)
gr.peakranger_dt_h3k4me3 <- sort(gr.peakranger_dt_h3k4me3)


# retain only unique regions, i.e. exclude common between DT and G1 samples
gr.peakranger_dt_h2az.unique <- gr.peakranger_dt_h2az[-which(gr.peakranger_dt_h2az$name %in% subsetByOverlaps(gr.peakranger_dt_h2az, gr.peakranger_g1_h2az)$name)]
gr.peakranger_g1_h2az.unique <- gr.peakranger_g1_h2az[-which(gr.peakranger_g1_h2az$name %in% subsetByOverlaps(gr.peakranger_g1_h2az, gr.peakranger_dt_h2az)$name)]
gr.peakranger_g1_h3k4me3.unique <- gr.peakranger_g1_h3k4me3[-which(gr.peakranger_g1_h3k4me3$name %in% subsetByOverlaps(gr.peakranger_g1_h3k4me3, gr.peakranger_dt_h3k4me3)$name)]
gr.peakranger_dt_h3k4me3.unique <- gr.peakranger_dt_h3k4me3[-which(gr.peakranger_dt_h3k4me3$name %in% subsetByOverlaps(gr.peakranger_dt_h3k4me3, gr.peakranger_g1_h3k4me3)$name)]

# results of using the ranger function (for narrow peaks)
dat.peakranger_narrow <- data.frame(wd = rep("/home/skurscheid/Data/Tremethick/LINE_1_project/peakranger_analysis/", 4), 
                                    cell = c("DT", "DT", "G1", "G1"), 
                                    mark = c("H2AZ", "H3K4me3", "H2AZ", "H3K4me3"), 
                                    peak = rep("narrow", 4))

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

#--------------------------------------------------------------------------
# comparison between MACS2 and PeakRanger
# G1 cells, H2A.Z
dat.macs2 <- data.frame(wd = as.character(rep("/home/skurscheid/Data/Tremethick/LINE_1_project/macs2_analysis/", 4)),
                        cell = c("DT", "DT", "G1", "G1"),
                        mark = c("H2AZ", "H3K4me3", "H2AZ", "H3K4me3"),
                        peak = rep("narrow", 4),
                        suffix = rep("_peaks.narrowPeak", 4))

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
grl.consensus <- lapply(seq_along(1:4), function(x) {
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
#   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19   X   Y   M 
# 144 173  70  79 109  82 112  54  56  43 172  34  84  39  37  42  79  39  42   7   0   0 
table(seqnames(gr.common.h2az.LINE1.flank1kb))
#  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  X  Y  M 
# 44 45  7  9 26 18 17 15 11 10 33  5 23  0  7  9 16  8  2  0  0  0 
# so peaks on 14 & X are missing
gr.consensus.LINE1.dt_h2az <- c(gr.consensus.LINE1.dt_h2az, grl.consensus.LINE1.flank1kb[[1]][which(seqnames(grl.consensus.LINE1.flank1kb[[1]]) %in% c("14", "X"))])
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
# so peaks on 14 & X are missing
gr.consensus.LINE1.g1_h2az <- c(gr.consensus.LINE1.g1_h2az, grl.consensus.LINE1.flank1kb[[3]][which(seqnames(grl.consensus.LINE1.flank1kb[[3]]) %in% c("14", "X"))])
table(seqnames(gr.consensus.LINE1.g1_h2az))
gr.consensus.LINE1.g1_h2az <- sort(gr.consensus.LINE1.g1_h2az)
grl.consensus.LINE1.g1_h2az <- GRangesList(lapply(canonicalChr, function(x) gr.consensus.LINE1.g1_h2az[which(seqnames(gr.consensus.LINE1.g1_h2az) == x)]))
irl.consensus.LINE1.g1_h2az <- as(grl.consensus.LINE1.g1_h2az, "IRangesList")

# import BAM files for visualization and quantitative analysis
# run picard MarkDuplicates prior to import to remove PCR duplicates
ga.dt_input <- import("DT_Input.sorted.rmdup.bam")
ga.dt_h2az <- import("DT_H2A.Z.sorted.rmdup.bam")
ga.g1_input <- import("G1_Input.sorted.rmdup.bam")
ga.g1_h2az <- import("G1_H2A.Z.sorted.rmdup.bam")

seqlevels(ga.dt_input, force = TRUE) <- canonicalChr
seqlevels(ga.dt_hs2az, force = TRUE) <- canonicalChr
seqlevels(ga.g1_input, force = TRUE) <- canonicalChr
seqlevels(ga.g1_h2az, force = TRUE) <- canonicalChr

Views(dt_h2az.cov, irl.consensus.LINE1.dt_h2az)


# regions unique to either cell line
gr.consensus_dt_h2az.unique <- gr.consensus_dt_h2az[-which(gr.consensus_dt_h2az$name %in% subsetByOverlaps(gr.consensus_dt_h2az, gr.consensus_g1_h2az)$name)]

