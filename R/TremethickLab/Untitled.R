#  TF1_H2A.Bbd_peak_analysis.R
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
library("TxDb.Hsapiens.UCSC.hg19.knownGene")
library("rtracklayer")
library("AnnotationDbi")


gr.macs.tf1 <- GRanges(seqnames = macs2.tf1$V1, IRanges(start = macs2.tf1$V2, end = macs2.tf1$V3), strand = "*", mcols = macs2.tf1[,c("V4", "V5", "V6", "V7", "V8", "V9", "V10")])
gr.macs.tf1a <- GRanges(seqnames = macs2.tf1a$V1, IRanges(start = macs2.tf1a$V2, end = macs2.tf1a$V3), strand = "*", mcols = macs2.tf1a[,c("V4", "V5", "V6", "V7", "V8", "V9", "V10")])

seqlevels(gr.macs.tf1, force = T) <- canonicalChr
seqlevels(gr.macs.tf1, force = T) <- canonicalChr
