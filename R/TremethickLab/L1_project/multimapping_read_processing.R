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

#----------class definitions--------------------------
RepeatRead <- setRefClass("RepeatRead",
                          fields = list(repClass = "numeric",
                                        repFamily = "numeric"),
                          methods = list(
                            unique = function(x) {
                              if (sum(repClass) == 1) {
                                return(TRUE)
                              } else {
                                return(FALSE)
                              }
                            },
                            rep.Class = function(x) {
                              names(repClass[which(repClass == 1)])
                            },
                            rep.Family = function(x) {
                              names(repFamily[which(repFamily == 1)])
                            }
                          ))

#----------variables & annotations--------------------
canonicalChr <- c(seq(1,19), "X", "Y", "M")

repeatMasker <- read.table("~/Data/Annotations/GRCm38_mm10/repeatMaskerMM10.txt", header = T, as.is = T, sep = "\t")
repeatMasker$genoName <- gsub("chr", "", repeatMasker$genoName)
gr.repeatMasker <- GRanges(seqnames = repeatMasker$genoName, 
                           IRanges(start = repeatMasker$genoStart, 
                                   end = repeatMasker$genoEnd), 
                           strand = repeatMasker$strand, 
                           repeatMasker[,c("repName", "repClass", "repFamily", "repStart", "repEnd", "repLeft", "swScore", "milliDiv", "milliDel", "milliIns")])
seqlevels(gr.repeatMasker, force = T) <- canonicalChr 

# convert GRanges object into GRangesList by repeat Class
grl.repeatMasker <- GRangesList(lapply(unique(mcols(gr.repeatMasker)[,"repClass"]), function(x) {gr.repeatMasker[which(mcols(gr.repeatMasker)[,"repClass"] == x)]}))
# convert GRangesList into IRangesList object for faster computation
irl.repeatMasker <- as(grl.repeatMasker, "IRangesList")

#----------load alignments from BAM file--------------
file1 <- ("../G1_Input_multimapped.sorted.sub.bam")

# this would load rather a lot of data, might be better to reduce it to only the relevant field
param1 <- ScanBamParam(what = c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar", "seq", "qual"), tag = c("AS", "XS"))
param2 <- ScanBamParam(what = c("qname", "flag", "rname", "strand", "pos", "mapq"), tag = c("AS", "XS"), which = )

# load BAM file
bam1 <- readGAlignments(file1, param = param2)
gal.bam1 <- readGAlignmentsList(file1, param = param2)
seqlevels(bam1, force = T) <- canonicalChr
seqlevels(gal.bam1, force = T) <- canonicalChr
names(gal.bam1) <- unique(unlist(lapply(gal.bam1, function(x) elementMetadata(x)$qname)))


# retrieve multimapping reads
w1 <- which(width(gal.bam1@partitioning) > 1)
gal.bam1 <- gal.bam1[w1]


# create IRangesList elements from each alignment per read
# irl.mm <- lapply(glb, function(x){
#   chr <- unique(as.character(elementMetadata(x)$rname))
#   #print(length(chr))
#   grl <- GRangesList(lapply(chr, function(y) {
#     x1 <- x[which(seqnames(x) == y)]
#     elementMetadata(x1) <- NULL
#     gr <- as(x1, "GRanges")
#   }))
#   names(grl) <- chr
#   irl <- as(grl, "IRangesList")
# })

# create vector to keep track of repeat Class
vrepClass <- rep(0, length(unique(mcols(gr.repeatMasker)[,"repClass"])))
names(vrepClass) <- unique(mcols(gr.repeatMasker)[,"repClass"])
vrepFamily <- rep(0, length(unique(mcols(gr.repeatMasker)[,"repFamily"])))
names(vrepFamily) <- unique(mcols(gr.repeatMasker)[,"repFamily"])

gal.bam1.ss <- gal.bam1[1:15]
# use snowfall for parallel execution
library(snowfall)
# initialize cluster
ptm <- proc.time()
sfInit(parallel = T, cpus = 8)
sfExport("gal.bam1.ss")
sfExport("gr.repeatMasker")
sfExport("vrepClass")
sfExport("vrepFamily")
sfExport(RepeatRead)
sfLibrary(GenomicRanges)
sfLibrary(Rsamtools)

# first step - determine which read maps to LINE1 elements only
# second step - isolate alignments of those reads
# third step - randomly choose one alignment
x1 <- sfLapply(gal.bam1.ss, function(x) {
  df1 <- mcols(subsetByOverlaps(gr.repeatMasker, x, ignore.strand = T))
  vrepClass[which(vrepClass > 0)] <- 0
  vrepFamily[which(vrepFamily > 0)] <- 0
  readRepClass <- unique(df1$repClass)
  readRepFamily <- unique(df1$repFamily)
  vrepClass[readRepClass] <- 1
  vrepFamily[readRepFamily] <- 1
  r1$repClass <- vrepClass
  r1$repFamily <- vrepFamily
  return(r1)
  })

proc.time() - ptm

# retain only list elements which are not null
x.select <- sfLapply(x1, function(x) !is.null(x))
x1 <- x1[unlist(x.select)]



# stopping the cluster
sfStop()

# fourth step - write reduced multimapper file back (should it include uniquely aligned reads?)

# trying the BamViews functionality 
# from: http://www.bioconductor.org/packages/release/bioc/vignettes/Rsamtools/inst/doc/Rsamtools-Overview.pdf
bamExperiment <- list(description = "LINE-1 elements views on DT/S/G1 H2AZ ChIP-Seq data", created = date())
bv <- BamViews(file1, bamRanges = gr.repeatMasker, bamExperiment = bamExperiment)

