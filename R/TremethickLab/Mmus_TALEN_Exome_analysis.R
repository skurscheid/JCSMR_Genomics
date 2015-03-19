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





