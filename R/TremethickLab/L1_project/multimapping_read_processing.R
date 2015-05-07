#  multimapping_read_processing.R
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

library("Rsamtools")
library("GenomicRanges")
library("GenomicFeatures")
library("GenomicAlignments")


#----------variables & annotations--------------------
canonicalChr <- c(seq(1,19), "X", "Y", "M")

repeatMasker <- read.table("~/Data/Annotations/GRCm38_mm10/repeatMaskerMM10.txt", header = T, as.is = T, sep = "\t")
repeatMasker$genoName <- gsub("chr", "", repeatMasker$genoName)
gr.repeatMasker <- GRanges(seqnames = repeatMasker$genoName, IRanges(start = repeatMasker$genoStart, end = repeatMasker$genoEnd), strand = repeatMasker$strand, mcols = repeatMasker[,c("repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft")])
seqlevels(gr.repeatMasker, force = T) <- canonicalChr 

#----------load alignments from BAM file--------------
file1 <- ("../G1_Input_multimapped.sorted.sub.bam")
param1 <- ScanBamParam(what = c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar", "seq", "qual"))
bam1 <- readGAlignments(file1, param = param1)
seqlevels(bam1, force = T) <- canonicalChr 

