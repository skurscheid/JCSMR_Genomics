library(Gviz)
library(Rsamtools)

# HOXA region as an example
chr <- "chr6"
from <- 52153727
to <- 52272467

# ctnnd2 region (from manuscript)
chr <- "chr15"
from <- 30430671
to <- 30781129
to.shorter <- 30481129

gr.mm <- GRanges(seqnames = chr, IRanges(start = from, end = to), strand = '*')

refGenes <- UcscTrack(genome = "mm9", chromosome = chr, track = "refGene", from = from, to = to, trackType = "GeneRegionTrack", rstarts = "exonStarts", rends = "exonEnds", gene = "name", symbol = "name2", transcript = "name", strand = "strand", fill = "#8282d2",stacking = "dense", name = "Other RefSeq")
save(refGenes, file = "/Volumes//gduserv.anu.edu.au/Data//Annotations/refGenes.rda")

load("../Annotations/refGenes.rda")
alTrack.Hippo.Lap1 <- AlignmentsTrack("Lap1_Hippo_adult7_NoIndex_L001_sorted.bam", isPaired = TRUE)
alTrack.Hippo.polyA <- AlignmentsTrack("./tophat_res/accepted_hits.bam", paired = FALSE)

pdf("temp.pdf", height = 6, width = 12)
plotTracks(trackList = c(refGenes, alTrack.Hippo.Lap1, alTrack.Hippo.polyA))
dev.off()
