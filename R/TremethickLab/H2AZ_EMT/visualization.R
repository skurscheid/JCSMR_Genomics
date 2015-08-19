# load libraries
library("Gviz")
library("rtracklayer")
library("biomaRt")
library("GenomicRanges")
library("GenomicAlignments")

# additional functions
# function to calculate binned averages from a coverage RleList and GRanges object
# adapted from http://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf
binnedAverage <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  means_list <- lapply(names(numvar),
                       function(seqname) {
                         views <- Views(numvar[[seqname]],
                                        bins_per_chrom[[seqname]])
                         viewMeans(views)
                       })
  new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}

binnedSum <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  means_list <- lapply(names(numvar),
                       function(seqname) {
                         views <- Views(numvar[[seqname]],
                                        bins_per_chrom[[seqname]])
                         viewSums(views)
                       })
  new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}

# set WD
setwd("~/Data/Tremethick/MDCK_sure_select/BAMs")


#----------connect to Ensembl biomaRt for annotation data----------------------
dog <- useMart("ensembl", dataset = "cfamiliaris_gene_ensembl")
filters <- listFilters(dog)

# # EMT markers
# # mesenchymal genes HGNC Symbols
# # FN1
# # ZEB1
# # TGFb1 -> TGFB1I1
# # SPARC
# # TWIST2
# mesenchymalMarkers <- c("FN1", "ZEB1", "TGFB1I1", "SPARC", "TWIST2")
# mesenchymalMarkers.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand","hgnc_symbol"), filters = "hgnc_symbol", values = mesenchymalMarkers, dog)
# gr.mesenchymalMarkers <- GRanges(mesenchymalMarkers.tab$chromosome_name, IRanges(mesenchymalMarkers.tab$start_position, mesenchymalMarkers.tab$end_position), strand = "*", hgnc_symbol = mesenchymalMarkers.tab$hgnc_symbol)
# 
# # epithelial markers
# # CDH1
# # SPP1
# # FGFBP1
# # MMP9
# # EPCAM
# epithelialMarkers <- c("CDH1", "SPP1", "FGFBP1", "MMP9", "EPCAM")
# epithelialMarkers.tab <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol"), filters = "hgnc_symbol", values = epithelialMarkers, dog)
# gr.epithelialMarkers <- GRanges(epithelialMarkers.tab$chromosome_name, IRanges(epithelialMarkers.tab$start_position, epithelialMarkers.tab$end_position), strand = "*", hgnc_symbol = epithelialMarkers.tab$hgnc_symbol)
# 
# gr.which <- c(gr.mesenchymalMarkers, gr.epithelialMarkers)
# 
# which.tab <- rbind(epithelialMarkers.tab, mesenchymalMarkers.tab)
# which.tab$type <- c(rep("epithelial", 5), rep("mesenchymal", 5))
# 

designTab <- data.frame(chromosome_name = as.character(c("chr5", "chr32", "chr3", "chr24", "chr10", "chr37", "chr2", "chr1", "chr4", "chr25")), 
                        start_position = c(80797492, 11334119, 64429408, 33254255, 49465757, 22436951, 15282325, 112608937, 57639351, 49145461),
                        end_position = c(80777492, 11374119, 64469409, 33294255, 49505757, 22476951, 15322325, 112648937, 57679351, 49185460),
                        hgnc_symbol = c("CDH1", "SPP1", "FBFBP1", "MMP9", "EPCAM", "FN1", "ZEB1", "TGFb1", "SPARC", "TWIST2"),
                        marker = c(rep("epithelial", 5), rep("mesenchymal", 5)))
designTab$chromosome_name <- as.character(designTab$chromosome_name)
#----------run MarkDuplicates on BAM files-------------------------------------


#----------H2AZ ChIP data------------------------------------------------------
# load ChIP coverage from BAM files
flag <- scanBamFlag(isProperPair = T, isPaired = T, isDuplicate = NA)
SBParam <- ScanBamParam(flag = flag, simpleCigar = F, what = scanBamWhat()) #, which = gr.which)

h2az.wt.rep1 <- readGAlignmentPairs("H2AZ-WT-rep1_S1_L001_sorted_MkDup.bam", param = SBParam)
h2az.wt.rep2 <- readGAlignmentPairs("H2AZ-WT-rep2_S2_L001_sorted_MkDup.bam", param = SBParam)
h2az.tgfb.rep1 <- readGAlignmentPairs("H2AZ-TGFb-rep1_S3_L001_sorted_MkDup.bam", param = SBParam)
h2az.tgfb.rep2 <- readGAlignmentPairs("H2AZ-TGFb-rep2_S4_L001_sorted_MkDup.bam", param = SBParam)

# correcting chromosome names to be compatible with UCSC (used for gene annotation in Gviz tracks)
seqlevels(h2az.wt.rep1, force = TRUE) <- paste("chr", seqlevels(h2az.wt.rep1), sep = "")
seqlevels(h2az.wt.rep2, force = TRUE) <- paste("chr", seqlevels(h2az.wt.rep2), sep = "")
seqlevels(h2az.tgfb.rep1, force = TRUE) <- paste("chr", seqlevels(h2az.tgfb.rep1), sep = "")
seqlevels(h2az.tgfb.rep2, force = TRUE) <- paste("chr", seqlevels(h2az.tgfb.rep2), sep = "")

# scaling factor for coverage across on-target regions [RPM]
h2az.wt.rep1.scale <- 1000000/length(h2az.wt.rep1)
h2az.wt.rep2.scale <- 1000000/length(h2az.wt.rep2)
h2az.tgfb.rep1.scale <- 1000000/length(h2az.tgfb.rep1)
h2az.tgfb.rep2.scale <- 1000000/length(h2az.tgfb.rep2)

# calculate coverage from both replicates
cov.h2az.wt <- (coverage(h2az.wt.rep1) * h2az.wt.rep1.scale + coverage(h2az.wt.rep2) * h2az.wt.rep2.scale) / 2
cov.h2az.tgfb <- (coverage(h2az.tgfb.rep1) * h2az.tgfb.rep1.scale + coverage(h2az.tgfb.rep2) * h2az.tgfb.rep2.scale) / 2

# TODO:
# calculate dyad sums, i.e. midpoints of fragments (Druliner et al. 2015 suggest 100bp windows in 10bp increments)


#----------Input data----------------------------------------------------------
input.wt.rep1 <- readGAlignmentPairs("Input-TGFb-rep1_S7_L001_sorted_MkDup.bam", param = SBParam)
input.wt.rep2 <- readGAlignmentPairs("Input-TGFb-rep2_S8_L001_sorted_MkDup.bam", param = SBParam)
input.tgfb.rep1 <- readGAlignmentPairs("Input-WT-rep1_S5_L001_sorted_MkDup.bam", param = SBParam)
input.tgfb.rep2 <- readGAlignmentPairs("Input-WT-rep2_S6_L001_sorted_MkDup.bam", param = SBParam)

# scaling factor for coverage for on-target reads [RPM]
input.wt.rep1.scale <- 1000000/length(input.wt.rep1)
input.wt.rep2.scale <- 1000000/length(input.wt.rep2)
input.tgfb.rep1.scale <- 1000000/length(input.tgfb.rep1)
input.tgfb.rep2.scale <- 1000000/length(input.tgfb.rep2)

# correcting chromosome names to be compatible with UCSC (used for gene annotation in Gviz tracks)
seqlevels(input.wt.rep1, force = TRUE) <- paste("chr", seqlevels(input.wt.rep1), sep = "")
seqlevels(input.wt.rep2, force = TRUE) <- paste("chr", seqlevels(input.wt.rep2), sep = "")
seqlevels(input.tgfb.rep1, force = TRUE) <- paste("chr", seqlevels(input.tgfb.rep1), sep = "")
seqlevels(input.tgfb.rep2, force = TRUE) <- paste("chr", seqlevels(input.tgfb.rep2), sep = "")

# calculate coverage from both replicates
cov.input.wt <- (coverage(input.wt.rep1) * input.wt.rep1.scale + coverage(input.wt.rep2) * input.wt.rep2.scale) / 2
cov.input.tgfb <- (coverage(input.tgfb.rep1) * input.tgfb.rep1.scale + coverage(input.tgfb.rep2) * input.tgfb.rep2.scale) / 2

#----------preparing data for plotting of coverage--------------
# creating single-bp level windows for visualization
gr.which.tiles <- tile(gr.which, width = 1) # width = bin size
gr.which.tiles <- unlist(gr.which.tiles)
seqlevels(gr.which.tiles, force = T) <- paste("chr", seqlevels(gr.which.tiles), sep = "")
seqlevels(gr.which.tiles, force = T) <- seqlevels(gr.which.tiles)[order(seqlevels(gr.which.tiles))]

cov.input.wt <- cov.input.wt[which(names(cov.input.wt) %in% seqlevels(gr.which.tiles))]
cov.input.wt <- cov.input.wt[names(cov.input.wt)[order(names(cov.input.wt))]]

cov.input.tgfb <- cov.input.tgfb[which(names(cov.input.tgfb) %in% seqlevels(gr.which.tiles))]
cov.input.tgfb <- cov.input.tgfb[names(cov.input.tgfb)[order(names(cov.input.tgfb))]]

cov.h2az.wt <- cov.h2az.wt[which(names(cov.h2az.wt) %in% seqlevels(gr.which.tiles))]
cov.h2az.wt <- cov.h2az.wt[names(cov.h2az.wt)[order(names(cov.h2az.wt))]]

cov.h2az.tgfb <- cov.h2az.tgfb[which(names(cov.h2az.tgfb) %in% seqlevels(gr.which.tiles))]
cov.h2az.tgfb <- cov.h2az.tgfb[names(cov.h2az.tgfb)[order(names(cov.h2az.tgfb))]]

# calculate binned averages from the 1bp-resolution coverage objects
bA.cov.input.wt <- binnedAverage(gr.which.tiles, cov.input.wt, "mean")
bA.cov.input.tgfb <- binnedAverage(gr.which.tiles, cov.input.tgfb, "mean")
bA.cov.h2az.wt <- binnedAverage(gr.which.tiles, cov.h2az.wt, "mean")
bA.cov.h2az.tgfb <- binnedAverage(gr.which.tiles, cov.h2az.tgfb, "mean")


#----------Using Gviz for visualization----------------------------------------
dT.cov.input.wt <- DataTrack(bA.cov.input.wt, type = "h", col = "darkgreen", name = "Input WT [rpm]")
dT.cov.input.tgfb <- DataTrack(bA.cov.input.tgfb, type = "h", col = "lightgreen", name = "Input TGFb [rpm]")
dT.cov.h2az.wt <- DataTrack(bA.cov.h2az.wt, type = "h", col = "darkred", name = "H2AZ WT [rpm]")
dT.cov.h2az.tgfb <- DataTrack(bA.cov.h2az.tgfb, type = "h", col = "red", name = "H2AZ TGFb [rpm]")

png(filename = "Plots.png", width = 1278, height = 720, units = "px")
pT <- apply(designTab[4,], 1, function(x){
  chromosome(dT.cov.input.wt) <- x["chromosome_name"]
  chromosome(dT.cov.input.tgfb) <- x["chromosome_name"]
  chromosome(dT.cov.h2az.wt) <- x["chromosome_name"]
  chromosome(dT.cov.h2az.tgfb) <- x["chromosome_name"]
  
  max.y <- max(max(values(dT.cov.input.wt)), max(values(dT.cov.input.tgfb)), max(values(dT.cov.h2az.wt)), max(values(dT.cov.h2az.tgfb)))
  displayPars(dT.cov.input.wt) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.input.tgfb) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.h2az.wt) <- list(ylim = c(0,max.y))
  displayPars(dT.cov.h2az.tgfb) <- list(ylim = c(0,max.y))

  biomTrack <- BiomartGeneRegionTrack(genome = "canFam3", chromosome = as.character(x["chromosome_name"]), start = as.integer(x["start_position"]), end = as.integer(x["end_position"]), name = paste(x["hgnc_symbol"], sep = ""))
  pT <- plotTracks(list(biomTrack, dT.cov.input.wt, dT.cov.h2az.wt, dT.cov.input.tgfb, dT.cov.h2az.tgfb), chromosome = x["chromosome_name"]) #from = as(x["start_position"], "integer"), to = as(x["end_position"], "integer")
  return(pT)
})
dev.off()


x <- which.tab[1,]
