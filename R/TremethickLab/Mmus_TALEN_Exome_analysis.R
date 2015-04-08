#  Mmus_TALEN_Exome_analysis.R
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
library("GenomicFeatures")
library("GenomicAlignments")
library("rtracklayer")
library("AnnotationDbi")

# some "global" variables
canonicalChr <- c(seq(1:19), "X", "Y", "M")


# get file names for deletions
fn <- list.files(pattern = "*deletions")

# load data & create GRanges objects
dels <- lapply(fn, function(x) data.frame(read.table(x, header = F, as.is = T, sep = "\t")))
grl.dels <- GRangesList(lapply(dels, function(x) GRanges(seqnames = x[,1], IRanges(start = x[,2], end = x[,3]), strand = "*", mcols = x[,4:5])))
seqlevels(grl.dels, force = T) <- canonicalChr
grl.dels <- GRangesList(lapply(grl.dels, function(x) sort(x)))

tab1 <- data.frame(samples = c("NM1_G2_23_index65", "NM4_G1_28_index9", "NM4_G2_18_index3", "NM4_G3_32_index61"), deletions = c("sample1", "sample2", "sample3", "sample4"))
data_dir = "/Volumes/MHS/researchdata/JCSMR/TremethickLab/Mmus_TALEN_Exomes/talen"

l.data <- lapply(1:length(tab1$samples), function(x) {
  p1 <- ScanBamParam(what = scanBamWhat(), which = grl.dels[[x]], flag = scanBamFlag(isDuplicate = F, isProperPair = T))
  b1 <- import(paste(data_dir, "/Sample_", tab1$samples[x], "/", tab1$samples[x], ".bam", sep = ""), param = p1)
})

l.data <- lapply(1:length(loc), function(x) {
  p1 <- ScanBamParam(what = scanBamWhat(), which = flank(loc[x], 250, both = TRUE), flag = scanBamFlag(isDuplicate = F, isProperPair = T))
  b1 <- GAlignmentsList(sapply(tab1$samples, function(y) import(paste(data_dir, "/Sample_", y, "/", y, ".bam", sep = ""), param = p1)))
})

l.cs.flank <- lapply(1:length(loc), function(x) {
  region1 <- flank(loc[x], 250, both = F)
  region2 <- GRanges(as.character(seqnames(loc[x])), IRanges(start = end(loc[x]), width = 250), strand = "*")
  af1 <- lapply(tab1$samples, function(y) alphabetFrequencyFromBam(BamFile(paste(data_dir, "/Sample_", y, "/", y, ".bam", sep = "")), param=region1, baseOnly=TRUE))
  af2 <- lapply(tab1$samples, function(y) alphabetFrequencyFromBam(BamFile(paste(data_dir, "/Sample_", y, "/", y, ".bam", sep = "")), param=region2, baseOnly=TRUE))
  af1 <- Reduce("+", af1)
  af2 <- Reduce("+", af2)
  
  cm1a1 <- t(af1[ , DNA_BASES])
  cm1a2 <- t(af2[ , DNA_BASES])
  cbind(cm_to_cs(cm1a1), cm_to_cs(cm1a2))
})



