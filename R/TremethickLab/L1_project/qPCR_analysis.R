#  qPCR_analysis.R
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


#---------load libraries----------------------------------
require(gdata)
require(GenomicRanges)
require(GenomicAlignments)
require(Rsamtools)

#--------global variables & settings----------------------
canonicalChr <- c(seq(1,19), "X", "Y", "M")

#---------LINE-1 5' UTR amplicon locations----------------

#---------First 5UTR primer pair--------------------------
qPCR.targets <- read.xls("Genome Location of LINE Primers - edited.xlsx", sheet = 1, as.is = T)
qPCR.targets <- data.frame(as.character(qPCR.targets[,1]))
qPCR.targets <- data.frame(as.character(qPCR.targets[grep("chr", qPCR.targets[,1]), ]))

qPCR.targets$chr <- rep(NA, nrow(qPCR.targets))
qPCR.targets$strand <- rep(NA, nrow(qPCR.targets))
qPCR.targets$start <- rep(0, nrow(qPCR.targets))
qPCR.targets$end <- rep(0, nrow(qPCR.targets))
qPCR.targets[,1] <- as.character(qPCR.targets[,1])
qPCR.targets$chr <- gsub("chr", "", unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[1])))
qPCR.targets$start <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[1])))
qPCR.targets$end <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[2])))
qPCR.targets[grep("\\+", qPCR.targets[,1]),]$strand <- "+"
qPCR.targets[grep("\\-", qPCR.targets[,1]),]$strand <- "-"

gr.qPCR.LINE1.5UTR1 <- GRanges(qPCR.targets$chr, IRanges(start = qPCR.targets$start, end = qPCR.targets$end), strand = qPCR.targets$strand)
seqlevels(gr.qPCR.LINE1.5UTR1, force = T) <- canonicalChr
gr.qPCR.LINE1.5UTR1$id <- "5UTR1"

#---------Second 5UTR primer pair-------------
qPCR.targets <- read.xls("Genome Location of LINE Primers - edited.xlsx", sheet = 2, as.is = T)
qPCR.targets <- data.frame(as.character(qPCR.targets[,1]))
qPCR.targets <- data.frame(as.character(qPCR.targets[grep("chr", qPCR.targets[,1]), ]))
qPCR.targets$chr <- rep(NA, nrow(qPCR.targets))
qPCR.targets$strand <- rep(NA, nrow(qPCR.targets))
qPCR.targets$start <- rep(0, nrow(qPCR.targets))
qPCR.targets$end <- rep(0, nrow(qPCR.targets))
qPCR.targets[,1] <- as.character(qPCR.targets[,1])
qPCR.targets$chr <- gsub("chr", "", unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[1])))
qPCR.targets$start <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[1])))
qPCR.targets$end <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[2])))
qPCR.targets[grep("\\+", qPCR.targets[,1]),]$strand <- "+"
qPCR.targets[grep("\\-", qPCR.targets[,1]),]$strand <- "-"

gr.qPCR.LINE1.5UTR2 <- GRanges(qPCR.targets$chr, IRanges(start = qPCR.targets$start, end = qPCR.targets$end), strand = qPCR.targets$strand)
seqlevels(gr.qPCR.LINE1.5UTR2, force = T) <- canonicalChr
gr.qPCR.LINE1.5UTR2$id <- "5UTR2"

#---------Third 5UTR primer pair-------------
qPCR.targets <- read.xls("Genome Location of LINE Primers - edited.xlsx", sheet = 3, as.is = T)
qPCR.targets <- data.frame(as.character(qPCR.targets[,1]))
qPCR.targets <- data.frame(as.character(qPCR.targets[grep("chr", qPCR.targets[,1]), ]))
qPCR.targets$chr <- rep(NA, nrow(qPCR.targets))
qPCR.targets$strand <- rep(NA, nrow(qPCR.targets))
qPCR.targets$start <- rep(0, nrow(qPCR.targets))
qPCR.targets$end <- rep(0, nrow(qPCR.targets))
qPCR.targets[,1] <- as.character(qPCR.targets[,1])
qPCR.targets$chr <- gsub("chr", "", unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[1])))
qPCR.targets$start <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[1])))
qPCR.targets$end <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[2])))
qPCR.targets[grep("\\+", qPCR.targets[,1]),]$strand <- "+"
qPCR.targets[grep("\\-", qPCR.targets[,1]),]$strand <- "-"

gr.qPCR.LINE1.5UTR3 <- GRanges(qPCR.targets$chr, IRanges(start = qPCR.targets$start, end = qPCR.targets$end), strand = qPCR.targets$strand)
seqlevels(gr.qPCR.LINE1.5UTR3, force = T) <- canonicalChr
gr.qPCR.LINE1.5UTR3$id <- "5UTR3"

gr.qPCR.LINE1.5UTR <- reduce(c(gr.qPCR.LINE1.5UTR1, gr.qPCR.LINE1.5UTR2, gr.qPCR.LINE1.5UTR3))
gr.qPCR.LINE1.5UTR$id <- "5UTR"

#---------L1ORF1A primer pair-----------
qPCR.targets <- read.xls("Genome Location of LINE Primers - edited.xlsx", sheet = 4, as.is = T)
qPCR.targets <- data.frame(as.character(qPCR.targets[,1]))
qPCR.targets <- data.frame(as.character(qPCR.targets[grep("chr", qPCR.targets[,1]), ]))
qPCR.targets$chr <- rep(NA, nrow(qPCR.targets))
qPCR.targets$strand <- rep(NA, nrow(qPCR.targets))
qPCR.targets$start <- rep(0, nrow(qPCR.targets))
qPCR.targets$end <- rep(0, nrow(qPCR.targets))
qPCR.targets[,1] <- as.character(qPCR.targets[,1])
qPCR.targets$chr <- gsub("chr", "", unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[1])))
qPCR.targets$start <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[1])))
qPCR.targets$end <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[2])))
qPCR.targets[grep("\\+", qPCR.targets[,1]),]$strand <- "+"
qPCR.targets[grep("\\-", qPCR.targets[,1]),]$strand <- "-"

gr.qPCR.LINE1.L1ORF1A <- GRanges(qPCR.targets$chr, IRanges(start = qPCR.targets$start, end = qPCR.targets$end), strand = qPCR.targets$strand)
seqlevels(gr.qPCR.LINE1.L1ORF1A, force = T) <- canonicalChr
gr.qPCR.LINE1.L1ORF1A$id <- "L1ORF1A"

#---------L1ORF1B primer pair--------------
qPCR.targets <- read.xls("Genome Location of LINE Primers - edited.xlsx", sheet = 5, as.is = T)
qPCR.targets <- data.frame(as.character(qPCR.targets[,1]))
qPCR.targets <- data.frame(as.character(qPCR.targets[grep("chr", qPCR.targets[,1]), ]))
qPCR.targets$chr <- rep(NA, nrow(qPCR.targets))
qPCR.targets$strand <- rep(NA, nrow(qPCR.targets))
qPCR.targets$start <- rep(0, nrow(qPCR.targets))
qPCR.targets$end <- rep(0, nrow(qPCR.targets))
qPCR.targets[,1] <- as.character(qPCR.targets[,1])
qPCR.targets$chr <- gsub("chr", "", unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[1])))
qPCR.targets$start <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[1])))
qPCR.targets$end <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[2])))
qPCR.targets[grep("\\+", qPCR.targets[,1]),]$strand <- "+"
qPCR.targets[grep("\\-", qPCR.targets[,1]),]$strand <- "-"

gr.qPCR.LINE1.L1ORF1B <- GRanges(qPCR.targets$chr, IRanges(start = qPCR.targets$start, end = qPCR.targets$end), strand = qPCR.targets$strand)
seqlevels(gr.qPCR.LINE1.L1ORF1B, force = T) <- canonicalChr
gr.qPCR.LINE1.L1ORF1B$id <- "L1ORF1B"

# L1ORF2A primer pair
qPCR.targets <- read.xls("Genome Location of LINE Primers - edited.xlsx", sheet = 6, as.is = T, skip = 1)
qPCR.targets <- qPCR.targets[,-c(1:8, 13)]
colnames(qPCR.targets) <- c("chr", "strand", "start", "end")

gr.qPCR.LINE1.L1ORF2A <- GRanges(qPCR.targets$chr, IRanges(start = qPCR.targets$start, end = qPCR.targets$end), strand = qPCR.targets$strand)
seqlevels(gr.qPCR.LINE1.L1ORF2A, force = T) <- canonicalChr
gr.qPCR.LINE1.L1ORF2A$id <- "L1ORF2A"

# L1ORF2B primer pair
qPCR.targets <- read.xls("Genome Location of LINE Primers - edited.xlsx", sheet = 7, as.is = T)
qPCR.targets <- data.frame(as.character(qPCR.targets[,1]))
qPCR.targets <- data.frame(as.character(qPCR.targets[grep("chr", qPCR.targets[,1]), ]))
qPCR.targets$chr <- rep(NA, nrow(qPCR.targets))
qPCR.targets$strand <- rep(NA, nrow(qPCR.targets))
qPCR.targets$start <- rep(0, nrow(qPCR.targets))
qPCR.targets$end <- rep(0, nrow(qPCR.targets))
qPCR.targets[,1] <- as.character(qPCR.targets[,1])
qPCR.targets$chr <- gsub("chr", "", unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[1])))
qPCR.targets$start <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[1])))
qPCR.targets$end <- as.integer(unlist(lapply(strsplit(unlist(lapply(strsplit(unlist(lapply(strsplit(qPCR.targets[,1], "\\ "), function(x) x[1])), ":"), function(x) x[2])), "[+-]"), function(x) x[2])))
qPCR.targets[grep("\\+", qPCR.targets[,1]),]$strand <- "+"
qPCR.targets[grep("\\-", qPCR.targets[,1]),]$strand <- "-"

gr.qPCR.LINE1.L1ORF2B <- GRanges(qPCR.targets$chr, IRanges(start = qPCR.targets$start, end = qPCR.targets$end), strand = qPCR.targets$strand)
seqlevels(gr.qPCR.LINE1.L1ORF2B, force = T) <- canonicalChr
gr.qPCR.LINE1.L1ORF2B$id <- "L1ORF2B"

gr.qPCR.LINE1.L1ORF1 <- reduce(c(gr.qPCR.LINE1.L1ORF1A, gr.qPCR.LINE1.L1ORF1B))
gr.qPCR.LINE1.L1ORF1$id <- "L1ORF1"

gr.qPCR.LINE1.L1ORF2 <- reduce(c(gr.qPCR.LINE1.L1ORF2A, gr.qPCR.LINE1.L1ORF2B))
gr.qPCR.LINE1.L1ORF2$id <- "L1ORF2"

# L1Md_A2 primer pair
qPCR.targets <- read.xls("Genome Location of LINE Primers - edited.xlsx", sheet = 8, as.is = T, skip = 1)
qPCR.targets <- qPCR.targets[,-c(1:8, 13)]
colnames(qPCR.targets) <- c("chr", "strand", "start", "end")

gr.qPCR.LINE1.L1Md_A2 <- GRanges(qPCR.targets$chr, IRanges(start = qPCR.targets$start, end = qPCR.targets$end), strand = qPCR.targets$strand)
seqlevels(gr.qPCR.LINE1.L1Md_A2, force = T) <- canonicalChr
gr.qPCR.LINE1.L1Md_A2$id <- "L1Md_A2"

# combine all into one
gr.qPCR.LINE1 <- c(gr.qPCR.LINE1.5UTR, gr.qPCR.LINE1.L1ORF1, gr.qPCR.LINE1.L1ORF2, gr.qPCR.LINE1.L1Md_A2)

# change seqlevels to UCSC encoding
gr.qPCR.LINE1.ucsc <- gr.qPCR.LINE1
seqlevels(gr.qPCR.LINE1.ucsc) <- paste("chr", as.character(seqlevels(gr.qPCR.LINE1)), sep = "")
names(gr.qPCR.LINE1.ucsc) <- paste(gr.qPCR.LINE1.ucsc$id, seqnames(gr.qPCR.LINE1.ucsc), start(gr.qPCR.LINE1.ucsc), end(gr.qPCR.LINE1.ucsc), sep = "_")
gs1 <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr.qPCR.LINE1.ucsc)
export(gs1, "qPCRTargets_reduced.fasta")

#----------load alignments from BAM file--------------
G1_H2AZ <- ("G1_H2AZ_multimapped.sorted.bam")
G1_Input <- ("G1_Input_multimapped.sorted.bam")
DT_H2AZ <- ("DT_H2AZ.multimapped.sorted.bam")
DT_Input <- ("DT_Input.multimapped.sorted.bam")
S_H2AZ <- ("S_H2AZ_multimapped.sorted.bam")
S_Input <- ("S_Input_multimapped.sorted.bam")


# this would load rather a lot of data, might be better to reduce it to only the relevant field
param1 <- ScanBamParam(what = c("qname", "flag", "rname", "strand", "pos", "mapq"), which = gr.qPCR.LINE1)

# load BAM file
bam1 <- readGAlignments(G1_H2AZ, param = param1)
seqlevels(bam1, force = T) <- canonicalChr
gal.bam1 <- readGAlignmentsList(G1_H2AZ, param = param1)
seqlevels(gal.bam1, force = T) <- canonicalChr
names(gal.bam1) <- unique(unlist(lapply(gal.bam1, function(x) elementMetadata(x)$qname)))
length(gal.bam1)


