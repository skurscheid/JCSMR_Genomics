library(biomaRt)
library(XML)
library(GenomicRanges)
library(GenomicAlignments)

#---------custom functions----------------------------------------------------------
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
{s <- substring(s, 2); if(strict) tolower(s) else s},
sep = "", collapse = " " )
sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

.summarize <- function(y){
  x <- grl1[[y]]
  lib <- libSizes[y]
  tss_up1000dn50_highExp <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.which[nMarkers_high]), x)))[,1]) / lib * 1000000
  tss_up1000dn50_lowExp <- sum(data.frame(union = assay(summarizeOverlaps(reduce(gr.which[nMarkers_low]), x)))[,1]) / lib * 1000000
  df1 <- data.frame(tss_up1000dn50_highExp, tss_up1000dn50_lowExp)
  return(df1)
}

source("~/Dropbox/Development/GeneralPurpose/R/heatmap.3.R")
source("~/Dropbox/Development/GeneralPurpose/R/plotCoverage.R")
source("~/Dropbox/Development/GeneralPurpose/R/plotCoverageStrands.R")

#---------global settings/variables-------------------------------------------------
canonicalChr <- c(seq(1,19), "X", "Y", "M")


#---------use ENSEMBL biomaRt for annotation data-----------------------------------
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# list available filters
filters <- listFilters(mouse)
# list available attributes
mmus.attribs <- listAttributes(mouse)
pages <- attributePages(mouse)
hsap.attribs <- listAttributes(human)


#---------scrape Wikipedia entry for neuronal marker names--------------------------
html1 <- readHTMLTable("http://en.wikipedia.org/wiki/Neuronal_lineage_marker#Neural_stem_cells_markers", header = T, trim = T, as.is = T)
neuro.markers <- data.frame(html1[[2]])
neuro.markers[,1] <- as.character(neuro.markers[,1])
neuro.markers[,2] <- as.character(neuro.markers[,2])

mmus.ensembl_genes <- getBM(attributes = c("ensembl_gene_id", "description", "chromosome_name", "start_position", "end_position", "strand"), mart = mouse)
hsap.eg1 <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "description"), mart = human)
hsap.eg2 <- getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"), mart = human)
hsap.ensembl_genes <- merge(hsap.eg1, hsap.eg2, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.y = T)
rm(hsap.eg1, hsap.eg2) 

hsap.ensembl_genes <- getBM(attributes = c("ensembl_gene_id", "description", "chromosome_name", "start_position", "end_position", "strand"), mart = human)
hsap.transcripts <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id"), mart = human)

#---------scan Ensembl genes for neuronal markers-----------------------------------
l2 <- lapply(neuro.markers[,2], function(x){
  markers <- unlist(strsplit(x, "; "))
  l1 <- lapply(markers, function(y){
        hgnc_symbol_rows <- c(grep(y, hsap.ensembl_genes$hgnc_symbol), grep(capwords(y, strict = T), hsap.ensembl_genes$hgnc_symbol))
        description_rows <- c(grep(y, hsap.ensembl_genes$description), grep(capwords(y, strict = T), hsap.ensembl_genes$description))
        r <- unique(c(hgnc_symbol_rows, description_rows))
        mmus_homol <- unique(hsap.ensembl_genes[r,]$mmusculus_homolog_ensembl_gene)
  })
  names(l1) <- markers
  return(l1)
})

names(l2) <- neuro.markers[,1]

mmus.neuronalMarkers <- lapply(l2, function(x) {
  x1 <- as.character(unlist(x))
  x1 <- x1[-which(x1 == "")]
})

mmus.neuronalMarkers.IDs <- unique(as.character(unlist(mmus.neuronalMarkers)))
mmus.neuronalMarkersTranscripts <- getBM(c("ensembl_gene_id", "description", "ensembl_transcript_id"), filters = "ensembl_gene_id", values = mmus.neuronalMarkers.IDs, mart = mouse)

mmus.Transcripts <- unique(rownames(quantDF))
mmus.Transcripts <- getBM(c("ensembl_gene_id", "description", "ensembl_transcript_id"), filters = "ensembl_transcript_id", values = mmus.Transcripts, mart = mouse)
mmus.Transcripts.IDs <- unique(mmus.Transcripts$ensembl_gene_id)

gr.neuronalMarkers <- GRanges(mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.neuronalMarkers.IDs,]$chromosome_name, IRanges(mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.neuronalMarkers.IDs,]$start_position, mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.neuronalMarkers.IDs,]$end_position), strand = mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.neuronalMarkers.IDs,]$strand, mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.neuronalMarkers.IDs,]$ensembl_gene_id, mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.neuronalMarkers.IDs,]$description)
gr.Transcripts <- GRanges(mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.Transcripts.IDs,]$chromosome_name, IRanges(mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.Transcripts.IDs,]$start_position, mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.Transcripts.IDs,]$end_position), strand = mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.Transcripts.IDs,]$strand, mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.Transcripts.IDs,]$ensembl_gene_id, mmus.ensembl_genes[mmus.ensembl_genes$ensembl_gene_id %in% mmus.Transcripts.IDs,]$description)

seqlevels(gr.neuronalMarkers, force = T) <- canonicalChr
seqlevels(gr.Transcripts, force = T) <- canonicalChr
colnames(elementMetadata(gr.Transcripts)) <- c("ensembl_gene_id", "description")
names(gr.Transcripts) <- elementMetadata(gr.Transcripts)$ensembl_gene_id

#---------load kallisto quantification results of Mmus hippocampus RNA-Seq----------
flist <- system("find /Volumes/gduserv/Data/Tremethick/Brain/mmus_hippocampus/Sample_Hippo_polyA_[123] -name abundance.txt", intern = T)

quantL <- lapply(flist, function(x){
  if (file.exists(x)){
  df <- read.table(x, header = T)
  }
})

dim(quantL[[2]])

quantDF.brain <- data.frame(matrix(ncol = 3, nrow = nrow(quantL[[1]])))
rownames(quantDF.brain) <- quantL[[1]][,1]
quantDF.brain[,1] <- quantL[[1]][,5]
quantDF.brain[,2] <- quantL[[2]][,5]
quantDF.brain[,3] <- quantL[[3]][,5]

colnames(quantDF.brain) <- c("mmus_hippocampus_RNA-Seq1", "mmus_hippocampus_RNA-Seq2", "mmus_hippocampus_RNA-Seq3")
quantDF.brain$mean <- apply(quantDF.brain, 1, mean)



flist <- system("find /Volumes/gduserv/Data/Tremethick/Testes/mmus_testes/RNA-Seq[123]* -name abundance.txt", intern = T)

quantL <- lapply(flist, function(x){
  if (file.exists(x)){
    df <- read.table(x, header = T)
  }
})

dim(quantL[[2]])

quantDF.testes <- data.frame(matrix(ncol = 3, nrow = nrow(quantL[[1]])))
rownames(quantDF.testes) <- quantL[[1]][,1]
quantDF.testes[,1] <- quantL[[1]][,5]
quantDF.testes[,2] <- quantL[[2]][,5]
quantDF.testes[,3] <- quantL[[3]][,5]
colnames(quantDF.testes) <- c("mmus_hippocampus_RNA-Seq1", "mmus_hippocampus_RNA-Seq2", "mmus_hippocampus_RNA-Seq3")
quantDF.testes$mean <- apply(quantDF.testes, 1, mean)

i1 <- intersect(rownames(quantDF.brain), mmus.neuronalMarkersTranscripts$ensembl_transcript_id)
quantDF.brain[i1,]
nMquantDF <- quantDF.brain[i1, ]

# select genes with > median(quantDF$mean) expression
nMarkers_high <- unique(mmus.neuronalMarkersTranscripts[which(mmus.neuronalMarkersTranscripts$ensembl_transcript_id%in% rownames(nMquantDF[nMquantDF$mean > 12,])),]$ensembl_gene_id)
nMarkers_low <- unique(mmus.neuronalMarkersTranscripts[which(mmus.neuronalMarkersTranscripts$ensembl_transcript_id%in% rownames(nMquantDF[nMquantDF$mean <= 12,])),]$ensembl_gene_id)
nMarkers_high <- intersect(nMarkers_high, names(gr.neuronalMarkers))
nMarkers_low <- intersect(nMarkers_low, names(gr.neuronalMarkers))

transcripts_high <- unique(mmus.Transcripts[which(mmus.Transcripts$ensembl_transcript_id %in% rownames(quantDF[quantDF$mean > 12,])),]$ensembl_gene_id)
transcripts_low <- unique(mmus.Transcripts[which(mmus.Transcripts$ensembl_transcript_id %in% rownames(quantDF[quantDF$mean <= 12,])),]$ensembl_gene_id)
transcripts_high <- intersect(transcripts_high, names(gr.Transcripts))
transcripts_low <- intersect(transcripts_low, names(gr.Transcripts))
transcripts_low <- transcripts_low[-which(transcripts_low %in% intersect(transcripts_high, transcripts_low))]

# flanking regions of genes
gr.neuronalMarkers.1500 <- flank(gr.neuronalMarkers, 1500, both = T)
gr.neuronalMarkers.1500 <- sort(gr.neuronalMarkers)

# TSS +50/-1500
gr.which <- promoters(gr.neuronalMarkers, upstream = 1500, down = 50)
gr.which <- promoters(gr.neuronalMarkers, upstream = 1000, down = 1000)
gr.which <- promoters(gr.neuronalMarkers, upstream = 500, down = 500)

gr.which <- promoters(gr.Transcripts, upstream = 1000, downstream = 1000)

#---------loading ChIP-Seq data--------------
fileList <- system("find ~/Data/Tremethick/Brain/mmus_hippocampus -name mm10_unique.sorted.bam", intern = T)

param1 <- ScanBamParam(what = c("qname", "flag", "rname", "strand", "pos", "mapq"), tag = c("AS", "XS"), which = gr.which)

# load files into list structure
BAMList <- lapply(fileList, function(x) {
  ga.bam1 <- readGAlignments(x, param = param1)
  return(ga.bam1)
})

BAMList <- GAlignmentsList(BAMList)
names(BAMList) <- c("H2A.Lap1_1", "H2A.Lap1_2", "H2A.Z_2", "H2A.Z_1", "Input_sonicated", "Input_MNase")

# remove duplicated reads
BAMList <- GAlignmentsList(lapply(BAMList, function(x) {
  ga1 <- x[!duplicated(elementMetadata(x)$pos)]
}))

# remove paired-end data (to simplify mix of single/paired-end libraries for now)
BAMList_temp <- GAlignmentsList(lapply(BAMList[c("H2A.Lap1_2", "H2A.Z_2", "Input_MNase", "Input_sonicated")], function(x) {
  ga1 <- x[which(elementMetadata(x)$flag %in% c(83,99))]
}))

BAMList[c("H2A.Lap1_2", "H2A.Z_2", "Input_MNase", "Input_sonicated")] <- BAMList_temp[c("H2A.Lap1_2", "H2A.Z_2", "Input_MNase", "Input_sonicated")]
rm(BAMList_temp)

# convert to GRanges and extend reads to 150bp fragments
grl1 <- GRangesList(lapply(BAMList, function(x) {
  gr1 <- granges(x)
 # width(gr1) <- 150
  seqlevels(gr1, force = T) <- canonicalChr
  return(gr1)
}))

cov.all <- lapply(grl1, function(x) coverage(x))

# entering only mapped fragments, i.e. mapped read pair counts as one fragment
libSizes <- c("H2A.Lap1_1" = 24799588, "H2A.Lap1_2" = 192775777, "H2A.Z_1" = 25824324, "H2A.Z_2" = 208673157, "Input_MNase" = 216581629, "Input_sonicated" = 258632248)

gr.which.high <- gr.which[nMarkers_high]
gr.which.low <- gr.which[nMarkers_low]

gr.which.high <- gr.which[transcripts_high]
gr.which.low <- gr.which[transcripts_low]

grl.which.high <- GRangesList(sapply(canonicalChr, function(x) gr.which.high[which(seqnames(gr.which.high) %in% x)]))
grl.which.low <- GRangesList(sapply(canonicalChr, function(x) gr.which.low[which(seqnames(gr.which.low) %in% x)]))
grl.which <- GRangesList(sapply(canonicalChr, function(x) gr.which[which(seqnames(gr.which) %in% x)]))

irl.which.high <- as(grl.which.high, "IRangesList")
irl.which.low <- as(grl.which.low, "IRangesList")
irl.which <- as(grl.which, "IRangesList")

cov.which.high <- lapply(cov.all, function(x) {v1 <- Views(x, irl.which.high); return(v1)})
cov.which.low <- lapply(cov.all, function(x) {v1 <- Views(x, irl.which.low); return(v1)})
cov.which <- lapply(cov.all, function(x) {v1 <- Views(x, irl.which); return(v1)})

#---------looking at GFAP specifically------------------
gr.gfap <- promoters(gr.gfap, 500, 500)
grl.gfap <- GRangesList(sapply(canonicalChr, function(x) gr.gfap[which(seqnames(gr.gfap) %in% x)]))
irl.gfap <- as(grl.gfap, "IRangesList")
cov.gfap <- lapply(cov.all, function(x) {v1 <- Views(x, irl.gfap); return(v1)})


# Input coverage for all 
l1 <- lapply(cov.which[["Input_MNase"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["Input_MNase"] * 1000000))
})
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
Input <- t(df1)
Input.mean <- apply(Input, 2, mean)

# Input coverage for all 
l1 <- lapply(cov.which[["Input_sonicated"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["Input_MNase"] * 1000000))
})
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
Input.sonicated <- t(df1)
Input.sonicated.mean <- apply(Input.sonicated, 2, mean)

#---------all counts in list()-------------------------------------------------------
lCovNeuroMarkersHigh <- lapply(seq_along(cov.which.high), function(w) {
  n1 <- names(cov.which.high)[[w]]
  print(n1)
  x <- cov.which.high[[n1]]
  l1 <- lapply(x, function(y) {
    lapply(y, function(z){
      if (mean(as.integer(z) > 0.6)) {
        as.numeric(as.integer(z) / libSizes[[n1]] * 1000000)
      }
    })
  })
  s1 <- unlist(lapply(l1, function(x1) length(x1) > 0))
  l1 <- l1[s1]
  ncols <- sum(unlist(lapply(l1, function(x2) {length(x2)})))
  print(ncols)
  df <- data.frame(matrix(unlist(l1), nrow = ncols, byrow = T))
})


s1 <- unlist(lapply(lCovNeuroMarkersHigh[["H2A.Lap1_2"]], function(x) length(x) > 0))



lCovNeuroMarkersLow <- lapply(seq_along(cov.which.low), function(w) {
  n1 <- names(cov.which.low)[[w]]
  print(n1)
  x <- cov.which.low[[n1]]
  lapply(x, function(y) {
    lapply(y, function(z){
      if (mean(as.integer(z) > 0.6)) {
        as.numeric(as.integer(z) / libSizes[[n1]] * 1000000)
      }
    })
  })
})

lCovTranscriptsHigh <- lapply(seq_along(cov.which.high), function(w) {
  n1 <- names(cov.which.high)[[w]]
  print(n1)
  x <- cov.which.high[[n1]]
  lapply(x, function(y) {
    lapply(y, function(z){
      if (mean(as.integer(z) > 0.6)) {
        as.numeric(as.integer(z) / libSizes[[n1]] * 1000000)
      }
    })
  })
})


lCovTranscriptsLow <- lapply(seq_along(cov.which.low), function(w) {
  n1 <- names(cov.which.low)[[w]]
  print(n1)
  x <- cov.which.low[[n1]]
  lapply(x, function(y) {
    lapply(y, function(z){
      if (mean(as.integer(z) > 0.6)) {
        as.numeric(as.integer(z) / libSizes[[n1]] * 1000000)
      }
    })
  })
})




l1 <- lapply(cov.which.high[["H2A.Lap1_1"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) {
    if(mean(as.integer(y) > 0.6))
    data.frame(as.integer(y) / libSizes["H2A.Lap1_1"] * 1000000)
    })
  })
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
h2alap1.high <- t(df1)

l1 <- lapply(cov.which.high[["H2A.Lap1_2"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["H2A.Lap1_2"] * 1000000))
  })
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
h2alap2.high <- t(df1)

l1 <- lapply(cov.which.high[["H2A.Z_1"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["H2A.Z_1"] * 1000000))
  })
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
h2az1.high <- t(df1)

l1 <- lapply(cov.which.high[["H2A.Z_2"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["H2A.Z_2"] * 1000000))
  })
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
h2az2.high <- t(df1)

l1 <- lapply(cov.which.high[["Input_MNase"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["Input_MNase"] * 1000000))
  })
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
Input.high <- t(df1)

Input.high.sum <- apply(Input.high, 2, sum)
h2az.high.sum <- apply(rbind(apply(h2az2.high, 2, sum), apply(h2az1.high, 2, sum)), 2, sum)
h2alap.high.sum <- apply(rbind(apply(h2alap1.high, 2, sum), apply(h2alap2.high, 2, sum)), 2, sum)

Input.high.mean <- apply(Input.high, 2, mean)
h2az.high.mean <- apply(rbind(apply(h2az2.high, 2, mean), apply(h2az1.high, 2, mean)), 2, mean)
h2alap.high.mean <- apply(rbind(apply(h2alap1.high, 2, mean), apply(h2alap2.high, 2, mean)), 2, mean)

# low expressed
l1 <- lapply(cov.which.low[["H2A.Lap1_1"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["H2A.Lap1_1"] * 1000000))
  })
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
h2alap1.low <- t(df1)

l1 <- lapply(cov.which.low[["H2A.Lap1_2"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["H2A.Lap1_2"] * 1000000))
  })
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
h2alap2.low <- t(df1)

l1 <- lapply(cov.which.low[["H2A.Z_1"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["H2A.Z_1"] * 1000000))
  })
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
h2az1.low <- t(df1)

l1 <- lapply(cov.which.low[["H2A.Z_2"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["H2A.Z_2"] * 1000000))
  })
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
h2az2.low <- t(df1)

l1 <- lapply(cov.which.low[["Input_MNase"]], function(x) if (length(x) > 1) {
  lapply(x, function(y) data.frame(as.integer(y) / libSizes["Input_MNase"] * 1000000))
  })
s1 <- unlist(lapply(l1, function(x) !is.null(x)))
df1 <- data.frame(l1[s1])
Input.low <- t(df1)

Input.low.sum <- apply(Input.low, 2, sum)
h2az.low.sum <- apply(rbind(apply(h2az2.low, 2, sum), apply(h2az1.low, 2, sum)), 2, sum)
h2alap.low.sum <- apply(rbind(apply(h2alap1.low, 2, sum), apply(h2alap2.low, 2, sum)), 2, sum)

Input.low.mean <- apply(Input.low, 2, mean)
h2az.low.mean <- apply(rbind(apply(h2az2.low, 2, mean), apply(h2az1.low, 2, mean)), 2, mean)
h2alap.low.mean <- apply(rbind(apply(h2alap1.low, 2, mean), apply(h2alap2.low, 2, mean)), 2, mean)

allTranscripts <- list(Input.sonicated.mean, Input.mean, Input.low.mean, Input.high.mean, h2alap.low.mean, h2alap.high.mean, h2az.low.mean, h2az.high.mean)
names(allTranscripts) <- c("Input.sonicated.mean", "Input.mean", "Input.low.mean", "Input.high.mean", "h2alap.low.mean", "h2alap.high.mean", "h2az.low.mean", "h2az.high.mean")

neuronalTranscripts <- list(Input.sonicated.mean, Input.mean, Input.low.mean, Input.high.mean, h2alap.low.mean, h2alap.high.mean, h2az.low.mean, h2az.high.mean)
names(neuronalTranscripts) <- c("Input.sonicated.mean", "Input.mean", "Input.low.mean", "Input.high.mean", "h2alap.low.mean", "h2alap.high.mean", "h2az.low.mean", "h2az.high.mean")

save(neuronalTranscripts,  file = "~/Data/Tremethick/Brain/neuronalTranscripts.rda")

l1 <- neuronalTranscripts
l1 <- allTranscripts
pdf("Metagene_plots_allTranscripts.pdf")

pdf("~/Data/Tremethick/Brain/Metagene_plots_neuroMarkers_H2ALap1.pdf")
pdf("~/Data/Tremethick/Brain/Metagene_plots_allTranscripts_H2ALap1.pdf")

plot(c(-500:499), l1[["h2alap.low.mean"]] / l1[["Input.sonicated.mean"]],
     type = "l",
     col = "white",
     ylim = c(min(c(min(l1[["h2alap.low.mean"]] / l1[["Input.sonicated.mean"]]), 
                    min(l1[["h2alap.high.mean"]] / l1[["Input.sonicated.mean"]]))),
              max(c(max(l1[["h2alap.low.mean"]] / l1[["Input.sonicated.mean"]]), 
                    max(l1[["h2alap.high.mean"]] / l1[["Input.sonicated.mean"]])))),
     main = "H2A.Lap1 at +/-500bp of TSS",
     bty = "l")
legend("bottomrigh", c("Transcripts, high", "Transcripts, low"), col = c("green", "red"), lty = 1, lwd = 3, bty = "n", cex = 0.8)
lines(lowess(l1[["h2alap.low.mean"]] / l1[["Input.sonicated.mean"]] ~ c(-500:499), f = 0.1),col="red", lwd = 4)
lines(lowess(l1[["h2alap.high.mean"]] / l1[["Input.sonicated.mean"]] ~ c(-500:499), f = 0.1),col="green", lwd = 4)
abline(v = 1, lwd = 4, lty = 2, col = "darkgrey")
dev.off()

pdf("~/Data/Tremethick/Brain/Metagene_plots_neuroMarkers_H2AZ.pdf")
pdf("~/Data/Tremethick/Brain/Metagene_plots_allTranscripts_H2AZ.pdf")

plot(c(-500:499), l1[["h2az.low.mean"]] / l1[["Input.mean"]],
     type = "l",
     col = "white",
     ylim = c(min(c(min(l1[["h2az.low.mean"]] / l1[["Input.mean"]]),
                    min(l1[["h2az.high.mean"]] / l1[["Input.mean"]]))), 
              max(c(max(l1[["h2az.low.mean"]] / l1[["Input.mean"]]),
                    max(l1[["h2az.high.mean"]] / l1[["Input.mean"]])))),
     main = "H2A.Z at +/- 500bp of TSS",
     bty = "l")
legend("topright", c("Transcripts, high", "Transcripts, low"), col = c("green", "red"), lty = 1, lwd = 3, bty = "n", cex = 0.8)
lines(lowess(l1[["h2az.low.mean"]] / l1[["Input.mean"]] ~ c(-500:499), f = 0.1),col="red", lwd = 4)
lines(lowess(l1[["h2az.high.mean"]] / l1[["Input.mean"]] ~ c(-500:499), f = 0.1),col="green", lwd = 4)
abline(v = 1, lwd = 4, lty = 2, col = "darkgrey")
dev.off()

boxplot(log2(nMquantDF[nMquantDF$mean > 12,]$mean + 1), log2(nMquantDF[nMquantDF$mean <= 12,]$mean + 1))
        
nMarkers_high_exp <- sapply(nMarkers_high, function(x) {
  t <- unique(mmus.neuronalMarkersTranscripts[mmus.neuronalMarkersTranscripts$ensembl_gene_id == x, "ensembl_transcript_id"])
  m <- mean(nMquantDF[t,]$mean)
})

nMarkers_input <- apply(Input, 1, mean)
nMarkers_input_low <- apply(Input.low, 1, mean)
nMarkers_input_high <- apply(Input.high, 1, mean)
nMarkers_high_h2az <- apply(cbind(apply(h2az2.high, 1, mean), apply(h2az1.high, 1, mean)), 1, mean)
nMarkers_high_h2alap1 <- apply(cbind(apply(h2alap1.high, 1, mean), apply(h2alap2.high, 1, mean)), 1, mean)

cor.test(nMarkers_high_exp[!is.na(nMarkers_high_exp)], sample(nMarkers_high_h2az, 79))
cor(nMarkers_high_exp[!is.na(nMarkers_high_exp)], sample(nMarkers_high_h2alap1, 79))

nMarkers_low_exp <- sapply(nMarkers_low, function(x) {
  t <- unique(mmus.neuronalMarkersTranscripts[mmus.neuronalMarkersTranscripts$ensembl_gene_id == x, "ensembl_transcript_id"])
  m <- mean(nMquantDF[t,]$mean)
})
nMarkers_low_h2az <- apply(cbind(apply(h2az2.low, 1, mean), apply(h2az1.low, 1, mean)), 1, mean)
nMarkers_low_h2alap1 <- apply(cbind(apply(h2alap1.low, 1, mean), apply(h2alap2.low, 1, mean)), 1, mean)

cor.test(nMarkers_low_exp[!is.na(nMarkers_low_exp)], sample(nMarkers_low_h2az, 172))
cor.test(nMarkers_low_exp[!is.na(nMarkers_low_exp)], sample(nMarkers_low_h2alap1, 172))

pdf("~/Data/Tremethick/Brain/Boxplot_H2AZ_Lap1_occupancy_neuroMarkers.pdf")
boxplot(log2(nMarkers_low_h2az / nMarkers_input_low + 1), log2(nMarkers_high_h2az / nMarkers_input_high + 1), names = c("Low", "High"), notch = T, main = "H2A.Z [log2 FC ChIP/Input]")
boxplot(log2(nMarkers_low_h2alap1 / nMarkers_input_low + 1),  log2(nMarkers_high_h2alap1 / nMarkers_input_high + 1), names = c("Low", "High"), notch = T, main = "H2A.Lap1 [log2 FC ChIP/Input]")
dev.off()

#---------combined plot of RPM only - H2A.Lap1 -----------

pdf("~/Data/Tremethick/Brain/H2ALap1_rpm_5000TSS500_V2.pdf", width = 12, height = 8)
ymin <- min(c(min(neuronalTranscripts[["h2alap.low.mean"]]), 
               min(neuronalTranscripts[["h2alap.high.mean"]]),
               min(allTranscripts[["h2alap.low.mean"]]),
               min(allTranscripts[["h2alap.high.mean"]])
               ))
               
ymax <- max(c(max(neuronalTranscripts[["h2alap.low.mean"]]), 
              max(neuronalTranscripts[["h2alap.high.mean"]]),
              max(allTranscripts[["h2alap.low.mean"]]),
              max(allTranscripts[["h2alap.high.mean"]])
              ))
ymin <- 0
ymax <- ymax + 0.005
par(mfrow = c(1,2))
l1 <- neuronalTranscripts
plot(c(-500:499), l1[["h2alap.low.mean"]],
     type = "l",
     col = "white",
     ylim = c(ymin, ymax),
     main = "H2A.Lap1",
     bty = "l",
     ylab = "[mean reads-per-million]",
     xlab = "Distance from TSS [bp]",
     cex.axis = 2,
     axes = F)
axis(side = 1, lwd = 3, at = c(-500, -250, 0, 250, 500))
axis(side = 2, lwd = 3)
legend("bottomright", c("Neuromarkers, high", "Neuromarkers, low"), col = c("green", "red"), lty = 1, lwd = 3, bty = "n", cex = 0.8)
lines(lowess(l1[["h2alap.low.mean"]] ~ c(-500:499), f = 0.1),col="red", lwd = 4)
lines(lowess(l1[["h2alap.high.mean"]] ~ c(-500:499), f = 0.1),col="green", lwd = 4)
abline(v = 1, lwd = 4, lty = 2, col = "darkgrey")

l1 <- allTranscripts
plot(c(-500:499), l1[["h2alap.low.mean"]],
     type = "l",
     col = "white",
     ylim = c(ymin, ymax),
     main = "",
     bty = "l",
     ylab = "",
     xlab = "",
     axes = F)
axis(side = 1, lwd = 3, at = c(-500, -250, 0, 250, 500))
#axis(side = 2, lwd = 3)
legend("bottomright", c("Transcripts, high", "Transcripts, low"), col = c("yellow", "blue"), lty = 1, lwd = 3, bty = "n", cex = 0.8)
lines(lowess(l1[["h2alap.low.mean"]] ~ c(-500:499), f = 0.1),col="blue", lwd = 4)
lines(lowess(l1[["h2alap.high.mean"]] ~ c(-500:499), f = 0.1),col="yellow", lwd = 4)
abline(v = 1, lwd = 4, lty = 2, col = "darkgrey")
dev.off()

#---------combined plot of RPM only - H2A.Z -----------

pdf("~/Data/Tremethick/Brain/H2AZ_rpm_5000TSS500_V2.pdf", width = 12, height = 8)
ymin <- min(c(min(neuronalTranscripts[["h2az.low.mean"]]), 
              min(neuronalTranscripts[["h2az.high.mean"]]),
              min(allTranscripts[["h2az.low.mean"]]),
              min(allTranscripts[["h2az.high.mean"]])
))

ymax <- max(c(max(neuronalTranscripts[["h2az.low.mean"]]), 
              max(neuronalTranscripts[["h2az.high.mean"]]),
              max(allTranscripts[["h2az.low.mean"]]),
              max(allTranscripts[["h2az.high.mean"]])
))
ymin <- 0
ymax <- ymax + 0.005
par(mfrow = c(1,2))
l1 <- neuronalTranscripts
plot(c(-500:499), l1[["h2az.low.mean"]],
     type = "l",
     col = "white",
     ylim = c(ymin, ymax),
     main = "H2A.Z",
     bty = "l",
     ylab = "[mean reads-per-million]",
     xlab = "Distance from TSS [bp]",
     cex.axis = 2,
     axes = F)
axis(side = 1, lwd = 3, at = c(-500, -250, 0, 250, 500))
axis(side = 2, lwd = 3)
legend("bottomright", c("Neuromarkers, high", "Neuromarkers, low"), col = c("green", "red"), lty = 1, lwd = 3, bty = "n", cex = 0.8)
lines(lowess(l1[["h2az.low.mean"]] ~ c(-500:499), f = 0.1),col="red", lwd = 4)
lines(lowess(l1[["h2az.high.mean"]] ~ c(-500:499), f = 0.1),col="green", lwd = 4)
abline(v = 1, lwd = 4, lty = 2, col = "darkgrey")

l1 <- allTranscripts
plot(c(-500:499), l1[["h2az.low.mean"]],
     type = "l",
     col = "white",
     ylim = c(ymin, ymax),
     main = "",
     bty = "l",
     ylab = "",
     xlab = "",
     axes = F)
axis(side = 1, lwd = 3, at = c(-500, -250, 0, 250, 500))
#axis(side = 2, lwd = 3)
legend("bottomright", c("Transcripts, high", "Transcripts, low"), col = c("yellow", "blue"), lty = 1, lwd = 3, bty = "n", cex = 0.8)
lines(lowess(l1[["h2az.low.mean"]] ~ c(-500:499), f = 0.1),col="blue", lwd = 4)
lines(lowess(l1[["h2az.high.mean"]] ~ c(-500:499), f = 0.1),col="yellow", lwd = 4)
abline(v = 1, lwd = 4, lty = 2, col = "darkgrey")
dev.off()

#---------combine all transcripts and neuromakers in one H2A.Z plot--------------
pdf("~/Data/Tremethick/Brain/H2AZ_rpm_5000TSS500_V3.pdf", width = 8, height = 8)
ymin <- min(c(min(neuronalTranscripts[["h2az.low.mean"]]), 
              min(neuronalTranscripts[["h2az.high.mean"]]),
              min(allTranscripts[["h2az.low.mean"]]),
              min(allTranscripts[["h2az.high.mean"]])
))

ymax <- max(c(max(neuronalTranscripts[["h2az.low.mean"]]), 
              max(neuronalTranscripts[["h2az.high.mean"]]),
              max(allTranscripts[["h2az.low.mean"]]),
              max(allTranscripts[["h2az.high.mean"]])
))
ymin <- 0
ymax <- ymax + 0.005
l1 <- neuronalTranscripts
plot(c(-500:499), l1[["h2az.low.mean"]],
     type = "l",
     col = "white",
     ylim = c(ymin, ymax),
     main = "H2A.Z",
     bty = "l",
     ylab = "[mean reads-per-million]",
     xlab = "Distance from TSS [bp]",
     cex.axis = 2,
     axes = F)
axis(side = 1, lwd = 3, at = c(-500, -250, 0, 250, 500))
axis(side = 2, lwd = 3)
legend("bottomright", c("Neuromarkers, high", "Neuromarkers, low", "Transcripts, high", "Transcripts, low"), col = c("green", "red", "yellow", "blue"), lty = 1, lwd = 3, bty = "n", cex = 0.8)
lines(lowess(l1[["h2az.low.mean"]] ~ c(-500:499), f = 0.1),col="red", lwd = 4)
lines(lowess(l1[["h2az.high.mean"]] ~ c(-500:499), f = 0.1),col="green", lwd = 4)
l1 <- allTranscripts
lines(lowess(l1[["h2az.low.mean"]] ~ c(-500:499), f = 0.1),col="blue", lwd = 4)
lines(lowess(l1[["h2az.high.mean"]] ~ c(-500:499), f = 0.1),col="yellow", lwd = 4)
abline(v = 0, lwd = 4, lty = 2, col = "darkgrey")
dev.off()

#---------combine all transcripts and neuromakers in one H2A.Lap1 plot--------------
pdf("~/Data/Tremethick/Brain/H2ALap1_rpm_5000TSS500_V3.pdf", width = 8, height = 8)
ymin <- min(c(min(neuronalTranscripts[["h2alap.low.mean"]]), 
              min(neuronalTranscripts[["h2alap.high.mean"]]),
              min(allTranscripts[["h2alap.low.mean"]]),
              min(allTranscripts[["h2alap.high.mean"]])
))

ymax <- max(c(max(neuronalTranscripts[["h2alap.low.mean"]]), 
              max(neuronalTranscripts[["h2alap.high.mean"]]),
              max(allTranscripts[["h2alap.low.mean"]]),
              max(allTranscripts[["h2alap.high.mean"]])
))
ymin <- 0
ymax <- ymax + 0.005
l1 <- neuronalTranscripts
plot(c(-500:499), l1[["h2alap.low.mean"]],
     type = "l",
     col = "white",
     ylim = c(ymin, ymax),
     main = "H2A.Lap1",
     bty = "l",
     ylab = "[mean reads-per-million]",
     xlab = "Distance from TSS [bp]",
     cex.axis = 2,
     axes = F)
axis(side = 1, lwd = 3, at = c(-500, -250, 0, 250, 500))
axis(side = 2, lwd = 3)
legend("bottomright", c("Neuromarkers, high", "Neuromarkers, low", "Transcripts, high", "Transcripts, low"), col = c("green", "red", "yellow", "blue"), lty = 1, lwd = 3, bty = "n", cex = 0.8)
lines(lowess(l1[["h2alap.low.mean"]] ~ c(-500:499), f = 0.1),col="red", lwd = 4)
lines(lowess(l1[["h2alap.high.mean"]] ~ c(-500:499), f = 0.1),col="green", lwd = 4)
l1 <- allTranscripts
lines(lowess(l1[["h2alap.low.mean"]] ~ c(-500:499), f = 0.1),col="blue", lwd = 4)
lines(lowess(l1[["h2alap.high.mean"]] ~ c(-500:499), f = 0.1),col="yellow", lwd = 4)
abline(v = 0, lwd = 4, lty = 2, col = "darkgrey")
dev.off()

#---------combine all transcripts and neuromakers in one Input MNase plot--------------
pdf("~/Data/Tremethick/Brain/Input_MNase_rpm_5000TSS500_V3.pdf", width = 8, height = 8)
ymin <- min(c(min(neuronalTranscripts[["Input.mean"]]), 
              min(allTranscripts[["Input.mean"]])
))

ymax <- max(c(max(neuronalTranscripts[["Input.mean"]]), 
              max(allTranscripts[["Input.mean"]])
))
ymin <- 0
l1 <- neuronalTranscripts
plot(c(-500:499), l1[["Input.mean"]],
     type = "l",
     col = "white",
     ylim = c(ymin, ymax),
     main = "Input, MNase",
     bty = "l",
     ylab = "[mean reads-per-million]",
     xlab = "Distance from TSS [bp]",
     cex.axis = 2,
     axes = F)
axis(side = 1, lwd = 3, at = c(-500, -250, 0, 250, 500))
axis(side = 2, lwd = 3)
legend("bottomright", c("Neuromarkers", "All Transcripts"), col = c("green", "yellow"), lty = 1, lwd = 3, bty = "n", cex = 0.8)
lines(lowess(l1[["Input.mean"]] ~ c(-500:499), f = 0.1),col="green", lwd = 4)
l1 <- allTranscripts
lines(lowess(l1[["Input.mean"]] ~ c(-500:499), f = 0.1),col="yellow", lwd = 4)
abline(v = 0, lwd = 4, lty = 2, col = "darkgrey")
dev.off()

#---------combine all transcripts and neuromakers in one Input MNase plot--------------
pdf("~/Data/Tremethick/Brain/Input_sonicated_rpm_5000TSS500_V3.pdf", width = 8, height = 8)
ymin <- min(c(min(neuronalTranscripts[["Input.sonicated.mean"]]), 
              min(allTranscripts[["Input.sonicated.mean"]])
))

ymax <- max(c(max(neuronalTranscripts[["Input.sonicated.mean"]]), 
              max(allTranscripts[["Input.sonicated.mean"]])
))
ymin <- 0
l1 <- neuronalTranscripts
plot(c(-500:499), l1[["Input.sonicated.mean"]],
     type = "l",
     col = "white",
     ylim = c(ymin, ymax),
     main = "Input, Sonicated",
     bty = "l",
     ylab = "[mean reads-per-million]",
     xlab = "Distance from TSS [bp]",
     cex.axis = 2,
     axes = F)
axis(side = 1, lwd = 3, at = c(-500, -250, 0, 250, 500))
axis(side = 2, lwd = 3)
legend("bottomright", c("Neuromarkers", "All Transcripts"), col = c("green", "yellow"), lty = 1, lwd = 3, bty = "n", cex = 0.8)
lines(lowess(l1[["Input.sonicated.mean"]] ~ c(-500:499), f = 0.1),col="green", lwd = 4)
l1 <- allTranscripts
lines(lowess(l1[["Input.sonicated.mean"]] ~ c(-500:499), f = 0.1),col="yellow", lwd = 4)
abline(v = 0, lwd = 4, lty = 2, col = "darkgrey")
dev.off()

#---------Loading MACS2 peak calls & converting into GRanges----------------
h2azPeaks <- read.table("~/Data/Tremethick/Brain/macs2_H2AZ/NA_peaks.xls", as.is = T, skip = 25, head = T)
h2alap1Peaks <- read.table("~/Data/Tremethick/Brain/macs2_H2Lap1/NA_peaks.xls", as.is = T, skip = 25, head = T)
gr.h2alap1Peaks <- GRanges(h2alap1Peaks$chr, IRanges(h2alap1Peaks$start, h2alap1Peaks$end), strand = "*", h2alap1Peaks$abs_summit, h2alap1Peaks$X.log10.pvalue., h2alap1Peaks$fold_enrichment, h2alap1Peaks$X.log10.qvalue.)
gr.h2azPeaks <- GRanges(h2azPeaks$chr, IRanges(h2azPeaks$start, h2azPeaks$end), strand = "*", h2azPeaks$abs_summit, h2azPeaks$X.log10.pvalue., h2azPeaks$fold_enrichment, h2azPeaks$X.log10.qvalue.)
seqlevels(gr.h2azPeaks, force = T) <- canonicalChr
seqlevels(gr.h2alap1Peaks, force =  T) <- canonicalChr

quantDF.testes$ensembl_transcript_id <- rownames(quantDF.testes)

m2 <- merge(mmus.Transcripts[grep("H2A", mmus.Transcripts$description), c("description", "ensembl_transcript_id")],
            quantDF.testes[mmus.Transcripts[grep("H2A", mmus.Transcripts$description),]$ensembl_transcript_id, ],
            by.x = "ensembl_transcript_id",
            by.y = "ensembl_transcript_id")
quantDF.testes$ensembl_transcript_id <- rownames(quantDF.testes)

quantDF.testes

SRR085725[hsap.transcripts[hsap.transcripts$ensembl_gene_id == "ENSG00000143013",]$ensembl_transcript_id,]

geneIDs <- c(LMO4 = "ENSG00000143013", #LMO4
             FOXP2 = "ENSG00000128573", #FOXP2
             H2Afb3 = "ENSG00000277745",
             H2Afb2 = "ENSG00000274183",
             H2Afb1 = "ENSG00000277858",
             GFAP = "ENSG00000131095", #GFAP
             H2AZ = "ENSG00000164032",
             H2AX = "ENSG00000188486",
             H2AY = "ENSG00000164032")

s1 <- unlist(sapply(geneIDs, function(x) hsap.transcripts[hsap.transcripts$ensembl_gene_id == x, ]$ensembl_transcript_id))

d1 <- SRR085474[s1,]
d1 <- rbind(d1, SRR085471[s1,])
d1 <- rbind(d1, SRR085725[s1,])
d1$gene <- rep(unlist(sapply(names(geneIDs), function(x) rep(x, length(grep(x, names(s1)))))),3)
#d1$gene <- as.factor(d1$gene)
#d1[d1$tpm == 0, ]$tpm <- NA

d2 <- summarySE(d1, measurevar = "tpm", groupvars = c("target_id", "gene"))
d2 <- d2[order(d2$gene),]

pdf("Human_brain_RNA-Seq.pdf", width = 12, height = 8)
p1 <- ggplot(d2, aes(x = target_id, y = tpm, fill = gene))
p1 <- p1 + geom_bar(stat="identity", colour = "grey")
p1 <- p1 + geom_errorbar(aes(ymin = tpm-se, ymax=tpm+se), size = .3, width = .2) 
p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
p1 <- p1 + facet_grid(. ~ gene, scales = "free", space = "free")
p1 <- p1 + theme(strip.text.x = element_text(angle = 90))
p1 <- p1 + ggtitle("Human Brain RNA-Seq data [three studies]")
p1
dev.off()

#--------------mouse---------------------
SRR287805 <- read.table("/Volumes/gduserv/Data/Tremethick/Brain/SRA_mmus_brain/SRR287805/abundance.txt", header = T, as.is = T, sep = "\t")
SRR287812 <- read.table("/Volumes/gduserv/Data/Tremethick/Brain/SRA_mmus_brain/SRR287812/abundance.txt", header = T, as.is = T, sep = "\t")
SRR287815 <- read.table("/Volumes/gduserv/Data/Tremethick/Brain/SRA_mmus_brain/SRR287815/abundance.txt", header = T, as.is = T, sep = "\t")
SRR287816 <- read.table("/Volumes/gduserv/Data/Tremethick/Brain/SRA_mmus_brain/SRR287816/abundance.txt", header = T, as.is = T, sep = "\t")

rownames(SRR287805) <- SRR287805$target_id
rownames(SRR287812) <- SRR287812$target_id
rownames(SRR287815) <- SRR287815$target_id
rownames(SRR287816) <- SRR287816$target_id

mmus.ensembl_genes[grep("foxp", mmus.ensembl_genes$description),]

mmus.geneIDs <- c(GFAP = "ENSMUSG00000020932",
                  FOXP2 = "ENSMUSG00000029563",
                  H2Afb1 = "ENSMUSG00000062651",
                  H2Afb2 = "ENSMUSG00000082482",
                  H2Afb3 = "ENSMUSG00000083616",
                  LMO4 = "ENSMUSG00000028266",
                  H2Afx = "ENSMUSG00000049932",
                  H2Afy = "ENSMUSG00000015937",
                  H2Afz = "ENSMUSG00000037894")

mmus.s1 <- unlist(sapply(mmus.geneIDs, function(x) mmus.transcripts[mmus.transcripts$ensembl_gene_id == x, ]$ensembl_transcript_id))

d1 <- SRR287805[mmus.s1,]
d1 <- rbind(d1, SRR287812[mmus.s1,])
d1 <- rbind(d1, SRR287815[mmus.s1,])
d1 <- rbind(d1, SRR287816[mmus.s1,])

d1$gene <- rep(unlist(sapply(names(mmus.geneIDs), function(x) rep(x, length(grep(x, names(mmus.s1)))))),4)
#d1$gene <- as.factor(d1$gene)
#d1[d1$tpm == 0, ]$tpm <- NA

d2 <- summarySE(d1, measurevar = "tpm", groupvars = c("target_id", "gene"))
d2 <- d2[order(d2$gene),]

pdf("Mouse_brain_RNA-Seq.pdf", width = 12, height = 8)
p1 <- ggplot(d2, aes(x = target_id, y = tpm, fill = gene))
p1 <- p1 + geom_bar(stat="identity", colour = "grey")
p1 <- p1 + geom_errorbar(aes(ymin = tpm-se, ymax=tpm+se), size = .3, width = .2) 
p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
p1 <- p1 + facet_grid(. ~ gene, scales = "free", space = "free")
p1 <- p1 + theme(strip.text.x = element_text(angle = 90))
p1 <- p1 + ggtitle("Mouse Brain RNA-Seq data [four studies]")
p1
dev.off()

