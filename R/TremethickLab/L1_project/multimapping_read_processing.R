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

library(Rsamtools)
library(GenomicRanges)
library(GenomicFeatures)
library(GenomicAlignments)


#----------variables & annotations--------------------
canonicalChr <- c(seq(1,19), "X", "Y", "M")

repeatMasker <- read.table("~/Data/Annotations/GRCm38_mm10/repeatMaskerMM10.txt", header = T, as.is = T, sep = "\t")
repeatMasker$genoName <- gsub("chr", "", repeatMasker$genoName)
gr.repeatMasker <- GRanges(seqnames = repeatMasker$genoName, IRanges(start = repeatMasker$genoStart, end = repeatMasker$genoEnd), strand = repeatMasker$strand, mcols = repeatMasker[,c("repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft")])
seqlevels(gr.repeatMasker, force = T) <- canonicalChr 

#----------load alignments from BAM file--------------
file1 <- ("../G1_Input_multimapped.sorted.sub.bam")

# this would load rather a lot of data, might be better to reduce it to only the relevant field
param1 <- ScanBamParam(what = c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar", "seq", "qual"), tag = c("AS", "XS"))
param2 <- ScanBamParam(what = c("qname", "flag", "rname", "strand", "pos", "mapq"), tag = c("AS", "XS"))

# load BAM files as GAlignmentsList object - every list entry is a read with one or more mappings
bam1 <- readGAlignmentsList(file1, param = param2)
seqlevels(bam1, force = T) <- canonicalChr

# retrieve multimapping reads
bam1.mm <- bam1[which(width(bam1@partitioning) > 1)]
#mm.id <- sapply(bam1.mm, function(x) unique(elementMetadata(x)$qname))

# use snowfall for parallel execution
library(snowfall)
# initialize cluster
ptm <- proc.time()
sfInit(parallel = T, cpus = 48)
sfExport("bam1.mm")
sfExport("gr.repeatMasker")
sfLibrary(GenomicRanges)
sfLibrary(Rsamtools)

# first step - determine which read maps to LINE1 elements only
# second step - isolate alignments of those reads
# third step - randomly choose one alignment
x1 <- sfLapply(bam1.mm, function(x) {
  t1 <- as.data.frame(mcols(subsetByOverlaps(gr.repeatMasker, x, ignore.strand = T)))[,c("mcols.repClass", "mcols.repFamily")]
  if (length(unique(t1$mcols.repClass)) == 1) {
    REclass <- unique(t1$mcols.repClass)
    if (REclass == "LINE") {
      x <- sample(x, 1)
      return(x)
    }
  }
})

# retain only list elements which are not null
x.select <- sfLapply(x1, function(x) !is.null(x))
x1 <- x1[unlist(x.select)]

proc.time() - ptm

# stopping the cluster
sfStop()

# fourth step - write reduced multimapper file back (should it include uniquely aligned reads?)


